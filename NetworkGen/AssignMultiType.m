function TypeData = AssignMultiType(obj, Atoms, Bonds)
% -------------------------------------------------------------------------
% AssignMultiType
%   Assign exported atom and bond types after CleanupNetwork.
%
%   This stage does not mutate Atoms or Bonds. It computes runtime labels
%   used by WriteDataFiles so the core generation pipeline can keep using
%   its existing internal matrix layout.
%
%   Connectivity rules are supplied as:
%       obj.architecture.types.connectivity = [type_i, type_j, bond_type, allowed]
%
%   Unspecified triples default to allowed. Atom types are assigned on the
%   cleaned graph using the atom-pair projection of those rules, and bond
%   types are then assigned to each surviving bond under the full
%   (atom-type, bond-type) connectivity table. When requested atom or bond
%   fractions cannot be achieved exactly while preserving the rules, the
%   closest valid assignment found is used and logged.
%
% INPUT
%   obj   : network object
%   Atoms : cleaned atom array
%   Bonds : cleaned bond array
%
% OUTPUT
%   TypeData : struct with atom_types, bond_types, and summary counts
% -------------------------------------------------------------------------

    TypeData = [];

    cfg = obj.architecture.types;
    if ~isfield(cfg, 'enabled') || ~logical(cfg.enabled)
        return;
    end

    natom = size(Atoms, 1);
    nbond = size(Bonds, 1);

    natom_type = max(1, round(cfg.natom_type));
    nbond_type = max(1, round(cfg.nbond_type));

    TypeData = struct(...
        'enabled',              true, ...
        'natom_type',           natom_type, ...
        'nbond_type',           nbond_type, ...
        'atom_types',           ones(natom, 1), ...
        'bond_types',           ones(nbond, 1), ...
        'atom_target_count',    zeros(1, natom_type), ...
        'atom_realized_count',  zeros(1, natom_type), ...
        'bond_target_count',    zeros(1, nbond_type), ...
        'bond_realized_count',  zeros(1, nbond_type), ...
        'atom_exact',           true, ...
        'bond_exact',           true, ...
        'connectivity_matrix',  true(natom_type, natom_type), ...
        'connectivity_tensor',  true(natom_type, natom_type, nbond_type) ...
    );

    obj.log.print('   [AssignMultiType] Assigning exported atom/bond types\n');

    atom_target = resolve_target_counts(natom, natom_type, cfg.atype_mode, ...
        cfg.atom_count, cfg.atom_frac);
    bond_target = resolve_target_counts(nbond, nbond_type, cfg.btype_mode, ...
        cfg.bond_count, cfg.bond_frac);

    TypeData.atom_target_count = atom_target;
    TypeData.bond_target_count = bond_target;

    [pair_allowed, bond_allowed] = resolve_connectivity_rules( ...
        cfg.connectivity, natom_type, nbond_type, obj);
    TypeData.connectivity_matrix = pair_allowed;
    TypeData.connectivity_tensor = bond_allowed;

    has_rule_restrictions = ~all(bond_allowed(:));

    if natom == 0
        obj.log.print('   [AssignMultiType] No atoms remain after cleanup; nothing to label\n');
        return;
    end

    if natom_type == 1
        atom_types = ones(natom, 1);
    elseif isempty(Bonds)
        atom_types = random_labels_from_counts(atom_target);
    elseif ~has_rule_restrictions
        atom_types = random_labels_from_counts(atom_target);
    else
        [atom_types, ~, ~] = assign_atom_types_constrained( ...
            Bonds, atom_target, bond_target, pair_allowed, bond_allowed, obj);
    end

    TypeData.atom_types = atom_types(:);
    TypeData.atom_realized_count = type_counts(TypeData.atom_types, natom_type);
    TypeData.atom_exact = isequal(TypeData.atom_realized_count, atom_target);

    if nbond == 0
        bond_types = zeros(0, 1);
        bond_realized = zeros(1, nbond_type);
    elseif nbond_type == 1
        bond_types = ones(nbond, 1);
        bond_realized = type_counts(bond_types, nbond_type);
    elseif ~has_rule_restrictions
        if isfield(cfg, 'btype_sel_method') && ~isempty(cfg.btype_sel_method) && ...
                ~strcmpi(cfg.btype_sel_method, 'random')
            obj.log.print(['   [AssignMultiType] Unsupported btype_sel_method="%s"; ' ...
                           'using random assignment\n'], cfg.btype_sel_method);
        end

        bond_types = random_labels_from_counts(bond_target);
        bond_realized = type_counts(bond_types, nbond_type);
    else
        [bond_types, bond_realized, n_unassigned_bonds, ~] = assign_bond_types_for_atoms( ...
            Bonds(:, 2), Bonds(:, 3), atom_types, bond_target, bond_allowed);

        if n_unassigned_bonds > 0
            error(['AssignMultiType: at least one cleaned bond has no allowed bond type ' ...
                   'under the requested connectivity rules.']);
        end
    end

    TypeData.bond_types = bond_types(:);
    TypeData.bond_realized_count = bond_realized;
    TypeData.bond_exact = isequal(TypeData.bond_realized_count, bond_target);

    obj.log.print('   [AssignMultiType] Atom target counts:  %s\n', mat2str(atom_target));
    obj.log.print('   [AssignMultiType] Atom realized counts:%s\n', mat2str(TypeData.atom_realized_count));

    if ~TypeData.atom_exact
        obj.log.print(['   [AssignMultiType] Exact atom targets were infeasible under ' ...
                       'the connectivity rules; using closest valid assignment\n']);
    end

    obj.log.print('   [AssignMultiType] Bond target counts:  %s\n', mat2str(bond_target));
    obj.log.print('   [AssignMultiType] Bond realized counts:%s\n', mat2str(TypeData.bond_realized_count));

    if ~TypeData.bond_exact
        obj.log.print(['   [AssignMultiType] Exact bond targets were infeasible under ' ...
                       'the connectivity rules; using closest valid assignment\n']);
    end

    record_type_stats(obj, TypeData);
end


function target_counts = resolve_target_counts(total_count, ntypes, mode, count_values, frac_values)

    if total_count <= 0
        target_counts = zeros(1, ntypes);
        return;
    end

    switch lower(mode)
        case {'fixed', 'count'}
            target_counts = reshape(count_values, 1, []);
            if numel(target_counts) ~= ntypes
                error('AssignMultiType: expected %d count entries, got %d.', ...
                    ntypes, numel(target_counts));
            end

            target_counts = max(0, round(target_counts));
            total_requested = sum(target_counts);

            if total_requested == total_count
                return;
            elseif total_requested == 0
                target_counts = distribute_by_largest_remainder( ...
                    (total_count / ntypes) * ones(1, ntypes), total_count);
            else
                scaled = (double(target_counts) / total_requested) * total_count;
                target_counts = distribute_by_largest_remainder(scaled, total_count);
            end

        case {'frac', 'fraction'}
            frac_values = reshape(frac_values, 1, []);
            if numel(frac_values) ~= ntypes
                error('AssignMultiType: expected %d fraction entries, got %d.', ...
                    ntypes, numel(frac_values));
            end

            frac_values = max(0, double(frac_values));
            frac_sum = sum(frac_values);
            if frac_sum <= 0
                error('AssignMultiType: fraction inputs must sum to a positive value.');
            end

            frac_values = frac_values / frac_sum;
            target_counts = distribute_by_largest_remainder(frac_values * total_count, total_count);

        otherwise
            error('AssignMultiType: unknown selection mode "%s".', mode);
    end
end


function counts = distribute_by_largest_remainder(real_values, total_count)

    real_values = reshape(real_values, 1, []);
    counts = floor(real_values);

    deficit = total_count - sum(counts);
    if deficit > 0
        remainders = real_values - counts;
        [~, order] = sort(remainders, 'descend');
        for kk = 1:deficit
            counts(order(kk)) = counts(order(kk)) + 1;
        end
    elseif deficit < 0
        remainders = real_values - counts;
        [~, order] = sort(remainders, 'ascend');
        for kk = 1:(-deficit)
            idx = order(kk);
            if counts(idx) > 0
                counts(idx) = counts(idx) - 1;
            end
        end
    end
end


function [pair_allowed, bond_allowed] = resolve_connectivity_rules(connectivity, natom_type, nbond_type, obj)

    bond_allowed = true(natom_type, natom_type, nbond_type);

    if isempty(connectivity)
        pair_allowed = any(bond_allowed, 3);
        return;
    end

    ncols = size(connectivity, 2);
    if ncols ~= 3 && ncols ~= 4
        error(['AssignMultiType: connectivity must be an N x 4 numeric array ' ...
               'of [type_i, type_j, bond_type, allowed]. ' ...
               'Legacy N x 3 [type_i, type_j, allowed] is also accepted.']);
    end

    if ncols == 3
        obj.log.print(['   [AssignMultiType] Interpreting legacy 3-column connectivity ' ...
                       'rules across all bond types\n']);
    end

    for row = 1:size(connectivity, 1)
        ti = round(connectivity(row, 1));
        tj = round(connectivity(row, 2));

        if ti < 1 || ti > natom_type || tj < 1 || tj > natom_type
            error('AssignMultiType: connectivity row %d references an out-of-range atom type.', row);
        end

        if ncols == 3
            tf = logical(round(connectivity(row, 3)));
            for bt = 1:nbond_type
                bond_allowed(ti, tj, bt) = tf;
                bond_allowed(tj, ti, bt) = tf;
            end
        else
            bt = round(connectivity(row, 3));
            tf = logical(round(connectivity(row, 4)));

            if bt < 1 || bt > nbond_type
                error('AssignMultiType: connectivity row %d references an out-of-range bond type.', row);
            end

            bond_allowed(ti, tj, bt) = tf;
            bond_allowed(tj, ti, bt) = tf;
        end
    end

    pair_allowed = any(bond_allowed, 3);
end


function labels = random_labels_from_counts(target_counts)

    total_count = sum(target_counts);
    labels = zeros(total_count, 1);
    ptr = 1;

    for tt = 1:numel(target_counts)
        nnow = target_counts(tt);
        if nnow <= 0
            continue;
        end

        labels(ptr:ptr+nnow-1) = tt;
        ptr = ptr + nnow;
    end

    if total_count > 1
        labels = labels(randperm(total_count));
    end
end


function [atom_types, bond_types, best_bond_realized] = assign_atom_types_constrained( ...
    Bonds, atom_target, bond_target, pair_allowed, bond_allowed, obj)

    natom = max(max(Bonds(:, 2:3), [], 'all'), 0);
    natom_type = numel(atom_target);

    if natom == 0
        atom_types = zeros(0, 1);
        bond_types = zeros(0, 1);
        best_bond_realized = zeros(1, numel(bond_target));
        return;
    end

    [edge_i, edge_j, neighbor_list, degree_vec] = build_graph(Bonds, natom);

    max_restarts = max(16, 6 * natom_type);
    best_assign = [];
    best_bond_types = [];
    best_bond_realized = zeros(1, numel(bond_target));
    best_conflicts = inf;
    best_unassigned = inf;
    best_bond_deviation = inf;
    best_deviation = inf;

    for attempt = 1:max_restarts
        [~, order] = sortrows([-degree_vec(:), rand(natom, 1)], [1 2]);
        assign = greedy_assignment(order, neighbor_list, atom_target, pair_allowed);
        assign = repair_rule_conflicts(assign, neighbor_list, atom_target, pair_allowed);

        [conflicts, realized] = assignment_score(assign, edge_i, edge_j, pair_allowed, natom_type);

        if conflicts == 0
            assign = rebalance_counts(assign, neighbor_list, atom_target, pair_allowed);
            [conflicts, realized] = assignment_score(assign, edge_i, edge_j, pair_allowed, natom_type);
        end

        [bond_types_now, bond_realized_now, n_unassigned, bond_deviation] = ...
            assign_bond_types_for_atoms(edge_i, edge_j, assign, bond_target, bond_allowed);

        deviation = sum(abs(realized - atom_target));

        if (conflicts < best_conflicts) || ...
           (conflicts == best_conflicts && n_unassigned < best_unassigned) || ...
           (conflicts == best_conflicts && n_unassigned == best_unassigned && ...
                bond_deviation < best_bond_deviation) || ...
           (conflicts == best_conflicts && n_unassigned == best_unassigned && ...
                bond_deviation == best_bond_deviation && deviation < best_deviation)
            best_assign = assign;
            best_bond_types = bond_types_now;
            best_bond_realized = bond_realized_now;
            best_conflicts = conflicts;
            best_unassigned = n_unassigned;
            best_bond_deviation = bond_deviation;
            best_deviation = deviation;
        end

        if best_conflicts == 0 && best_unassigned == 0 && ...
                best_bond_deviation == 0 && best_deviation == 0
            break;
        end
    end

    if isempty(best_assign) || best_conflicts > 0
        error(['AssignMultiType: could not find an atom-type assignment that ' ...
               'satisfies the requested connectivity rules on the cleaned network.']);
    end

    if best_unassigned > 0
        error(['AssignMultiType: could not find a valid atom-type assignment with at least ' ...
               'one allowed bond type for every cleaned bond.']);
    end

    atom_types = best_assign(:);
    bond_types = best_bond_types(:);
    if best_deviation > 0
        obj.log.print(['   [AssignMultiType] Connectivity rules constrained the ' ...
                       'requested atom fractions/counts\n']);
    end
end


function [bond_types, realized, n_unassigned, deviation] = assign_bond_types_for_atoms( ...
    edge_i, edge_j, atom_types, target_counts, bond_allowed)

    nbond = numel(edge_i);
    nbond_type = numel(target_counts);

    if nbond == 0
        bond_types = zeros(0, 1);
        realized = zeros(1, nbond_type);
        n_unassigned = 0;
        deviation = 0;
        return;
    end

    atom_i = atom_types(edge_i);
    atom_j = atom_types(edge_j);

    feasible = false(nbond, nbond_type);
    for bt = 1:nbond_type
        bt_vec = bt * ones(nbond, 1);
        lin_idx = sub2ind(size(bond_allowed), atom_i, atom_j, bt_vec);
        feasible(:, bt) = bond_allowed(lin_idx);
    end

    bond_types = zeros(nbond, 1);
    n_unassigned = sum(~any(feasible, 2));
    if n_unassigned > 0
        realized = zeros(1, nbond_type);
        deviation = inf;
        return;
    end

    eligible_counts = sum(feasible, 1);
    alt_counts = sum(feasible, 2);
    current_counts = zeros(1, nbond_type);
    remaining_target = target_counts;
    unassigned = true(nbond, 1);

    scarcity = eligible_counts - target_counts;
    [~, type_order] = sort(scarcity + 1e-6 * rand(1, nbond_type), 'ascend');

    for bt = type_order
        candidates = find(unassigned & feasible(:, bt));
        if isempty(candidates)
            continue;
        end

        need = max(0, remaining_target(bt));
        if need <= 0
            continue;
        end

        cand_score = alt_counts(candidates) + 1e-3 * rand(numel(candidates), 1);
        [~, order] = sort(cand_score, 'ascend');
        chosen = candidates(order(1:min(need, numel(candidates))));

        bond_types(chosen) = bt;
        unassigned(chosen) = false;
        current_counts(bt) = current_counts(bt) + numel(chosen);
        remaining_target(bt) = remaining_target(bt) - numel(chosen);
    end

    remaining_bonds = find(unassigned);
    rem_score = alt_counts(remaining_bonds) + 1e-3 * rand(numel(remaining_bonds), 1);
    [~, rem_order] = sort(rem_score, 'ascend');
    remaining_bonds = remaining_bonds(rem_order);

    for kk = 1:numel(remaining_bonds)
        bond_idx = remaining_bonds(kk);
        feasible_types = find(feasible(bond_idx, :));

        deficits = remaining_target(feasible_types);
        if any(deficits > 0)
            pick_score = deficits + 1e-6 * rand(1, numel(feasible_types));
            [~, idx] = max(pick_score);
        else
            overfill = current_counts(feasible_types) - target_counts(feasible_types);
            pick_score = overfill + 1e-6 * rand(1, numel(feasible_types));
            [~, idx] = min(pick_score);
        end

        bt = feasible_types(idx);
        bond_types(bond_idx) = bt;
        current_counts(bt) = current_counts(bt) + 1;
        remaining_target(bt) = remaining_target(bt) - 1;
        unassigned(bond_idx) = false;
    end

    [bond_types, current_counts] = rebalance_bond_counts( ...
        bond_types, feasible, target_counts, current_counts);

    realized = current_counts;
    deviation = sum(abs(realized - target_counts));
end


function [bond_types, counts_now] = rebalance_bond_counts(bond_types, feasible, target_counts, counts_now)

    max_passes = 10;

    for pass = 1:max_passes
        deficits = target_counts - counts_now;
        surplus = counts_now - target_counts;

        if ~any(deficits > 0) || ~any(surplus > 0)
            break;
        end

        changed = false;
        deficit_types = find(deficits > 0);
        [~, deficit_order] = sort(deficits(deficit_types), 'descend');
        deficit_types = deficit_types(deficit_order);

        for bt_def = deficit_types
            candidates = find((bond_types ~= bt_def) & feasible(:, bt_def));
            if isempty(candidates)
                continue;
            end

            current_type = bond_types(candidates);
            movable = surplus(current_type) > 0;
            candidates = candidates(movable);
            if isempty(candidates)
                continue;
            end

            cand_score = sum(feasible(candidates, :), 2) + 1e-3 * rand(numel(candidates), 1);
            [~, order] = sort(cand_score, 'ascend');
            chosen = candidates(order(1));

            bt_old = bond_types(chosen);
            bond_types(chosen) = bt_def;
            counts_now(bt_old) = counts_now(bt_old) - 1;
            counts_now(bt_def) = counts_now(bt_def) + 1;
            changed = true;
        end

        if ~changed
            break;
        end
    end
end


function assign = greedy_assignment(order, neighbor_list, target_counts, allowed)

    natom = numel(order);
    natom_type = numel(target_counts);

    assign = zeros(natom, 1);
    current_counts = zeros(1, natom_type);

    for kk = 1:natom
        node = order(kk);
        neigh = neighbor_list{node};
        neigh = neigh(assign(neigh) > 0);
        neigh_types = assign(neigh);

        scores = inf(1, natom_type);
        for tt = 1:natom_type
            local_conflicts = count_local_conflicts(tt, neigh_types, allowed);
            overfill = max(0, current_counts(tt) + 1 - target_counts(tt));
            shortage = max(0, target_counts(tt) - current_counts(tt));
            scores(tt) = 1e6 * local_conflicts + 1e3 * overfill - shortage + rand() * 1e-3;
        end

        [~, best_type] = min(scores);
        assign(node) = best_type;
        current_counts(best_type) = current_counts(best_type) + 1;
    end
end


function assign = repair_rule_conflicts(assign, neighbor_list, target_counts, allowed)

    natom_type = numel(target_counts);
    max_passes = 8;

    for pass = 1:max_passes
        changed = false;
        for node = randperm(numel(assign))
            cur_type = assign(node);
            neigh_types = assign(neighbor_list{node});
            cur_conflicts = count_local_conflicts(cur_type, neigh_types, allowed);

            compatible = fully_compatible_types(neigh_types, allowed, natom_type);
            if ~any(compatible)
                continue;
            end

            counts_now = type_counts(assign, natom_type);
            best_type = cur_type;
            best_score = 1e3 * cur_conflicts + sum(abs(counts_now - target_counts));

            for tt = find(compatible).'
                trial_counts = counts_now;
                trial_counts(cur_type) = trial_counts(cur_type) - 1;
                trial_counts(tt) = trial_counts(tt) + 1;
                trial_score = sum(abs(trial_counts - target_counts));

                if cur_conflicts > 0
                    trial_score = trial_score - 1e2;
                end

                if (tt ~= cur_type) && (trial_score < best_score)
                    best_type = tt;
                    best_score = trial_score;
                end
            end

            if best_type ~= cur_type
                assign(node) = best_type;
                changed = true;
            end
        end

        if ~changed
            break;
        end
    end
end


function assign = rebalance_counts(assign, neighbor_list, target_counts, allowed)

    natom_type = numel(target_counts);
    max_passes = 10;

    for pass = 1:max_passes
        counts_now = type_counts(assign, natom_type);
        deficit_types = find(counts_now < target_counts);
        if isempty(deficit_types)
            break;
        end

        surplus_nodes = find(counts_now(assign).' > target_counts(assign).');
        if isempty(surplus_nodes)
            break;
        end

        moved = false;
        surplus_nodes = surplus_nodes(randperm(numel(surplus_nodes)));

        for node = surplus_nodes(:).'
            cur_type = assign(node);
            neigh_types = assign(neighbor_list{node});
            compatible = fully_compatible_types(neigh_types, allowed, natom_type);

            candidate_types = deficit_types(compatible(deficit_types));
            if isempty(candidate_types)
                continue;
            end

            deficits = target_counts(candidate_types) - counts_now(candidate_types);
            [~, idx] = max(deficits);
            new_type = candidate_types(idx);

            if new_type ~= cur_type
                assign(node) = new_type;
                moved = true;
                break;
            end
        end

        if ~moved
            break;
        end
    end
end


function [conflicts, realized] = assignment_score(assign, edge_i, edge_j, allowed, natom_type)

    if isempty(edge_i)
        conflicts = 0;
    else
        lin_idx = sub2ind(size(allowed), assign(edge_i), assign(edge_j));
        conflicts = sum(~allowed(lin_idx));
    end

    realized = type_counts(assign, natom_type);
end


function conflicts = count_local_conflicts(type_id, neigh_types, allowed)

    if isempty(neigh_types)
        conflicts = 0;
        return;
    end

    conflicts = sum(~allowed(type_id, neigh_types));
end


function compatible = fully_compatible_types(neigh_types, allowed, natom_type)

    if isempty(neigh_types)
        compatible = true(natom_type, 1);
        return;
    end

    compatible = all(allowed(:, neigh_types), 2);
end


function counts = type_counts(labels, ntypes)

    if isempty(labels)
        counts = zeros(1, ntypes);
        return;
    end

    labels = labels(:);
    labels = labels(labels >= 1 & labels <= ntypes);

    if isempty(labels)
        counts = zeros(1, ntypes);
        return;
    end

    counts = accumarray(double(labels), 1, [ntypes, 1], @sum, 0).';
end


function [edge_i, edge_j, neighbor_list, degree_vec] = build_graph(Bonds, natom)

    edge_i = Bonds(:, 2);
    edge_j = Bonds(:, 3);
    neighbor_list = cell(natom, 1);

    for kk = 1:size(Bonds, 1)
        ii = edge_i(kk);
        jj = edge_j(kk);
        neighbor_list{ii}(end+1) = jj; %#ok<AGROW>
        neighbor_list{jj}(end+1) = ii; %#ok<AGROW>
    end

    degree_vec = zeros(natom, 1);
    if ~isempty(edge_i)
        degree_vec = degree_vec + accumarray(double([edge_i; edge_j]), 1, [natom, 1], @sum, 0);
    end
end


function record_type_stats(obj, TypeData)

    obj.log.record('multitype_enabled', true, ...
        'desc', 'Post-cleanup multi-type export enabled');

    for tt = 1:numel(TypeData.atom_realized_count)
        obj.log.record(sprintf('atom_type_%d_count', tt), TypeData.atom_realized_count(tt), ...
            'desc', sprintf('Atom count for type %d', tt));
    end

    for tt = 1:numel(TypeData.bond_realized_count)
        obj.log.record(sprintf('bond_type_%d_count', tt), TypeData.bond_realized_count(tt), ...
            'desc', sprintf('Bond count for type %d', tt));
    end

    obj.log.record('multitype_atom_targets_exact', TypeData.atom_exact, ...
        'desc', 'Exact atom type targets achieved');
    obj.log.record('multitype_bond_targets_exact', TypeData.bond_exact, ...
        'desc', 'Exact bond type targets achieved');
end