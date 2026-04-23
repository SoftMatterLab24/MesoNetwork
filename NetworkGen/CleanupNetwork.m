function [Atoms, Bonds, Nvec] = CleanupNetwork(obj, Atoms, Bonds, Nvec)
% -------------------------------------------------------------------------
% CleanupNetwork
%   Post-processing cleanup after bond generation and defect insertion.
%
%   Performs the following steps in order:
%     1. Prune isolated atoms (degree = 0)
%     2. Iteratively prune low-degree atoms below obj.peratom.min_degree_keep
%     3. Identify all connected components via union-find
%     4. WARN if multiple large components exist
%     5. Discard all components except the largest
%     6. Renumber atom IDs and bond endpoints consecutively from 1
%     7. Rebuild per-atom neighbor lists
%     8. Filter Nvec to match surviving bonds (if non-empty)
%
% Atoms column layout (new):
%   [ID | molID | X | Y | Z | deg | nbr_1 ... nbr_{Max_peratom_bond}]
%
% INPUT
%   obj   : network object
%   Atoms : atom array
%   Bonds : [M x 5] bond array [id | i | j | L0 | type]
%   Nvec  : per-bond quantity array, or [] when not yet assigned
%
% OUTPUT
%   Atoms, Bonds, Nvec : cleaned and renumbered
% -------------------------------------------------------------------------

    LARGE_COMP_BOND_THRESHOLD = 20;

    if isempty(Atoms) || isempty(Bonds)
        obj.log.print('   [CleanupNetwork] Empty Atoms or Bonds on entry; nothing to clean.\n');
        Atoms = zeros(0, max(6 + obj.peratom.Max_peratom_bond, 6));
        Bonds = zeros(0, 5);
        Nvec  = [];
        return;
    end

    min_keep         = obj.peratom.min_degree_keep;
    Max_peratom_bond = obj.peratom.Max_peratom_bond;

    % For bottle-brush geometry, disable minimum degree pruning: rod atoms
    % (especially rod ends) have low degree from just rod bonds.
    if strcmpi(obj.architecture.geometry, 'bottle_brush')
        min_keep = 0;
    end

    natom_in  = size(Atoms, 1);
    nbond_in  = size(Bonds, 1);

    obj.log.print('   [CleanupNetwork] Starting: %d atoms, %d bonds\n', natom_in, nbond_in);

    [Atoms, Bonds, Nvec, n_pruned_atoms, n_pruned_bonds] = ...
        iterative_degree_prune(Atoms, Bonds, Nvec, min_keep);

    if n_pruned_atoms > 0
        obj.log.print('   [CleanupNetwork] Degree pruning removed %d atoms and %d bonds\n', ...
            n_pruned_atoms, n_pruned_bonds);
    else
        obj.log.print('   [CleanupNetwork] Degree pruning: nothing removed\n');
    end

    if isempty(Atoms) || isempty(Bonds)
        warning('CleanupNetwork: all atoms/bonds removed during degree pruning.');
        Atoms = zeros(0, 6 + Max_peratom_bond);
        Bonds = zeros(0, 5);
        Nvec  = [];
        return;
    end

    %% ------------------------------------------------------------------
    %  Union-find to label connected components
    %% ------------------------------------------------------------------
    surv_ids = Atoms(:, 1);
    n_surv   = numel(surv_ids);

    id_to_loc = build_id_map(surv_ids);

    uf_parent = int32(1:n_surv);
    uf_rank   = zeros(n_surv, 1, 'int32');

    for k = 1:size(Bonds, 1)
        loc1 = id_to_loc(int64(Bonds(k, 2)));
        loc2 = id_to_loc(int64(Bonds(k, 3)));

        r1 = uf_find(loc1, uf_parent);
        r2 = uf_find(loc2, uf_parent);

        if r1 ~= r2
            [uf_parent, uf_rank] = uf_union(r1, r2, uf_parent, uf_rank);
        end
    end

    roots = zeros(n_surv, 1, 'int32');
    for ii = 1:n_surv
        roots(ii) = uf_find(ii, uf_parent);
    end

    unique_roots = unique(roots);
    n_comp       = numel(unique_roots);

    comp_atom_count = zeros(n_comp, 1);
    comp_bond_count = zeros(n_comp, 1);

    root_to_comp = containers.Map('KeyType','int32','ValueType','int32');
    for c = 1:n_comp
        root_to_comp(unique_roots(c)) = int32(c);
        comp_atom_count(c) = sum(roots == unique_roots(c));
    end

    for k = 1:size(Bonds, 1)
        loc1 = id_to_loc(int64(Bonds(k, 2)));
        c    = root_to_comp(roots(loc1));
        comp_bond_count(c) = comp_bond_count(c) + 1;
    end

    large_comp_mask  = comp_bond_count >= LARGE_COMP_BOND_THRESHOLD;
    n_large_comp     = sum(large_comp_mask);

    if n_large_comp > 1
        warning(['CleanupNetwork: %d disconnected regions each have >= %d bonds. ' ...
                 'The network is not well connected -- consider increasing atom ' ...
                 'density (rho_atom), reducing void coverage, or widening the ' ...
                 'bond cutoff radius.  Only the largest component will be kept.\n' ...
                 '   Component bond counts: %s'], ...
                 n_large_comp, LARGE_COMP_BOND_THRESHOLD, ...
                 num2str(sort(comp_bond_count(large_comp_mask), 'descend')'));
    end

    [~, main_comp_idx] = max(comp_atom_count);
    main_root          = unique_roots(main_comp_idx);
    in_main            = (roots == main_root);

    n_small_atoms = sum(~in_main);
    if n_small_atoms > 0
        Atoms     = Atoms(in_main, :);
        keep_ids  = surv_ids(in_main);

        keep_b = ismember(Bonds(:, 2), keep_ids) & ismember(Bonds(:, 3), keep_ids);
        n_small_bonds = sum(~keep_b);
        Bonds = Bonds(keep_b, :);

        if ~isempty(Nvec)
            Nvec = Nvec(keep_b, :);
        end

        obj.log.print(['   [CleanupNetwork] Removed %d atoms and %d bonds ' ...
                 'in %d small/secondary disconnected component(s)\n'], ...
            n_small_atoms, n_small_bonds, n_comp - 1);
    else
        obj.log.print('   [CleanupNetwork] Single connected component, no fragments removed\n');
    end

    if isempty(Atoms) || isempty(Bonds)
        warning('CleanupNetwork: no atoms/bonds remain after component pruning.');
        Atoms = zeros(0, 6 + Max_peratom_bond);
        Bonds = zeros(0, 5);
        Nvec  = [];
        return;
    end

    %% ------------------------------------------------------------------
    %  Renumber atom IDs and bond endpoints consecutively
    %% ------------------------------------------------------------------
    natom_new     = size(Atoms, 1);
    old_ids       = Atoms(:, 1);
    new_ids       = (1:natom_new)';
    Atoms(:, 1)   = new_ids;

    if ~isempty(Bonds)
        [~, loc1]    = ismember(Bonds(:, 2), old_ids);
        [~, loc2]    = ismember(Bonds(:, 3), old_ids);
        Bonds(:, 2)  = new_ids(loc1);
        Bonds(:, 3)  = new_ids(loc2);
        Bonds(:, 1)  = (1:size(Bonds, 1))';
    end

    %% ------------------------------------------------------------------
    %  Rebuild per-atom neighbor lists (degree at col 6, nbrs at 7..)
    %% ------------------------------------------------------------------
    needed_cols = 6 + Max_peratom_bond;
    cur_cols    = size(Atoms, 2);
    if cur_cols < needed_cols
        Atoms(:, cur_cols+1:needed_cols) = 0;
    end

    Atoms(:, 6)                       = 0;   % reset degree
    Atoms(:, 7:6 + Max_peratom_bond)  = 0;   % clear neighbor slots

    for k = 1:size(Bonds, 1)
        ii = Bonds(k, 2);
        jj = Bonds(k, 3);

        nb_i = Atoms(ii, 6) + 1;
        if nb_i <= Max_peratom_bond
            Atoms(ii, 6)        = nb_i;
            Atoms(ii, 6 + nb_i) = jj;
        end

        nb_j = Atoms(jj, 6) + 1;
        if nb_j <= Max_peratom_bond
            Atoms(jj, 6)        = nb_j;
            Atoms(jj, 6 + nb_j) = ii;
        end
    end

    obj.log.print('   [CleanupNetwork] Done.  %d atoms (-%d), %d bonds (-%d)\n', ...
        size(Atoms, 1), natom_in - size(Atoms, 1), ...
        size(Bonds, 1), nbond_in - size(Bonds, 1));

end


%% =======================================================================
%  LOCAL HELPER FUNCTIONS
%% =======================================================================

function [Atoms, Bonds, Nvec, n_pruned_atoms, n_pruned_bonds] = ...
        iterative_degree_prune(Atoms, Bonds, Nvec, min_keep)

    n_pruned_atoms = 0;
    n_pruned_bonds = 0;

    if isempty(Atoms) || isempty(Bonds)
        return;
    end

    thresh = max(1, min_keep);   % at least 1: removes isolated atoms
    if min_keep == 0
        % Allow degree-0 atoms to be kept (e.g. bottle_brush with rod bonds
        % only). Still remove truly isolated atoms implicitly by the
        % component step later on - here we just skip degree pruning.
        return;
    end

    changed = true;
    while changed
        natom = size(Atoms, 1);
        ids   = Atoms(:, 1);

        deg = zeros(natom, 1);
        if ~isempty(Bonds)
            id_map = build_id_map(ids);
            for k = 1:size(Bonds, 1)
                r1 = id_map(int64(Bonds(k, 2)));
                r2 = id_map(int64(Bonds(k, 3)));
                deg(r1) = deg(r1) + 1;
                deg(r2) = deg(r2) + 1;
            end
        end

        low_deg = deg < thresh;
        if ~any(low_deg)
            changed = false;
            break;
        end

        if ~isempty(Bonds)
            kill_ids = ids(low_deg);
            kill_b   = ismember(Bonds(:, 2), kill_ids) | ...
                       ismember(Bonds(:, 3), kill_ids);

            n_pruned_bonds = n_pruned_bonds + sum(kill_b);
            Bonds = Bonds(~kill_b, :);

            if ~isempty(Nvec)
                Nvec = Nvec(~kill_b, :);
            end
        end

        n_pruned_atoms = n_pruned_atoms + sum(low_deg);
        Atoms  = Atoms(~low_deg, :);
        changed = true;

        if isempty(Atoms) || isempty(Bonds)
            break;
        end
    end
end


function id_map = build_id_map(ids)
    n      = numel(ids);
    id_map = containers.Map('KeyType','int64','ValueType','int32');
    for ii = 1:n
        id_map(int64(ids(ii))) = int32(ii);
    end
end


function root = uf_find(x, uf_parent)
    while uf_parent(x) ~= x
        uf_parent(x) = uf_parent(uf_parent(x));
        x            = uf_parent(x);
    end
    root = x;
end


function [uf_parent, uf_rank] = uf_union(r1, r2, uf_parent, uf_rank)
    if uf_rank(r1) < uf_rank(r2)
        uf_parent(r1) = r2;
    elseif uf_rank(r1) > uf_rank(r2)
        uf_parent(r2) = r1;
    else
        uf_parent(r2)   = r1;
        uf_rank(r1)     = uf_rank(r1) + 1;
    end
end
