function [Atoms, Bonds] = AddBondsMono(obj, Atoms, LatticeData, PreBonds)
% -------------------------------------------------------------------------
% AddBondsMono
%
% Mono bond-topology generation:
%   - random      -> distance-limited random bonding
%   - hex_lattice -> strict first-neighbor lattice bonding
%
% Atoms column layout (read-only here except for degree/neighbor rebuild):
%   [ID | molID | X | Y | Z | deg | nbr_1 ... nbr_{Max_peratom_bond}]
%
% INPUT:
%   obj         : network object
%   Atoms       : atom array
%   LatticeData : lattice metadata (required for hex_lattice)
%   PreBonds    : optional pre-existing bond list to seed adjacency with so
%                 duplicate bonds aren't created. Only adjacency is seeded
%                 (NOT degree counters), so e.g. a mid-rod atom can still
%                 accept Max_peratom_bond crosslinks on top of its 2 rod
%                 bonds. When PreBonds is non-empty, internal pruning is
%                 skipped (CleanupNetwork handles it later).
%
% OUTPUT:
%   Atoms : updated atom array
%   Bonds : [bondID | atomID_i | atomID_j | L0 | type=0]
%           type column is left 0; AssignMultiType assigns it later.
% -------------------------------------------------------------------------

    if nargin < 4
        PreBonds = [];
    end

    geom = lower(obj.architecture.geometry);

    switch geom

        case 'random'
            [Atoms, Bonds] = connect_random_mono(obj, Atoms, PreBonds);

        case 'hex_lattice'
            [Atoms, Bonds] = connect_lattice_mono(obj, Atoms, LatticeData);

        case 'bottle_brush'
            % For bottle-brush we fall through to the random connector
            % and use PreBonds (rod bonds) to seed adjacency.
            [Atoms, Bonds] = connect_random_mono(obj, Atoms, PreBonds);

        otherwise
            error('AddBondsMono: unknown geometry "%s".', obj.architecture.geometry);
    end

end


% =========================================================================
% RANDOM MONO
% =========================================================================
function [AtomsOut, BondsOut] = connect_random_mono(obj, Atoms, PreBonds)

    natom             = size(Atoms,1);
    Max_bond          = obj.domain.Max_bond;
    Max_peratom_bond  = obj.peratom.Max_peratom_bond;
    global_limit      = obj.domain.bond_global_try_limit;
    stall_limit       = obj.domain.max_attempts_without_progress;
    min_keep          = obj.peratom.min_degree_keep;

    xlo = obj.domain.xlo; xhi = obj.domain.xhi;
    ylo = obj.domain.ylo; yhi = obj.domain.yhi;
    Lx = xhi - xlo; Ly = yhi - ylo;

    isPeriodic = strcmpi(obj.domain.boundary, 'periodic');

    a = obj.domain.min_node_sep;
    if strcmpi(obj.architecture.spacing_multiplier_mode, 'auto')
        spacing_multiplier = 1.8;
    else
        spacing_multiplier = obj.architecture.spacing_multiplier;
    end

    Rcut = a * spacing_multiplier;
    Rcut2 = Rcut * Rcut;

    ids = Atoms(:,1);
    x   = Atoms(:,3);   % X at col 3 under new layout
    y   = Atoms(:,4);   % Y at col 4

    hx = Rcut;
    hy = Rcut;
    nx = max(1, floor((xhi - xlo)/hx));
    ny = max(1, floor((yhi - ylo)/hy));

    cx = floor((x - xlo)/hx) + 1;
    cy = floor((y - ylo)/hy) + 1;
    cx = max(1, min(nx, cx));
    cy = max(1, min(ny, cy));

    Cells = cell(nx, ny);
    for i = 1:natom
        Cells{cx(i), cy(i)}(end+1) = i; %#ok<AGROW>
    end

    deg = zeros(natom,1);
    adj = sparse(natom, natom);

    % Seed adjacency from PreBonds so duplicate bonds aren't generated.
    % Degree is NOT seeded intentionally: atoms already carrying rod bonds
    % can still receive Max_peratom_bond crosslinks. Atom IDs in PreBonds
    % equal row indices at this point in the pipeline (no pruning yet).
    if ~isempty(PreBonds)
        for k = 1:size(PreBonds, 1)
            r1 = PreBonds(k, 2);
            r2 = PreBonds(k, 3);
            adj(r1, r2) = 1;
            adj(r2, r1) = 1;
        end
    end

    BondsRows = zeros(Max_bond, 3);   % [row1 row2 L]
    nbond = 0;
    ntries = 0;
    no_progress = 0;

    tic
    while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)

        ntries = ntries + 1;

        r1 = randi(natom);
        if deg(r1) >= Max_peratom_bond
            no_progress = no_progress + 1;
            continue;
        end

        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
        if isempty(neigh)
            no_progress = no_progress + 1;
            continue;
        end

        x1 = x(r1);
        y1 = y(r1);

        cand = zeros(16,1);
        ncan = 0;

        for kk = 1:numel(neigh)
            r2 = neigh(kk);

            if r2 == r1
                continue;
            end
            if deg(r2) >= Max_peratom_bond
                continue;
            end
            if adj(r1,r2) ~= 0
                continue;
            end

            dxv = x(r2) - x1;
            dyv = y(r2) - y1;
            d = minimum_image(isPeriodic, dxv, dyv, Lx, Ly);

            if (d*d) < Rcut2
                ncan = ncan + 1;
                if ncan > numel(cand)
                    cand = [cand; zeros(numel(cand),1)]; %#ok<AGROW>
                end
                cand(ncan) = r2;
            end
        end

        if ncan == 0
            no_progress = no_progress + 1;
            continue;
        end

        r2 = cand(randi(ncan));
        L = minimum_image(isPeriodic, x(r2)-x(r1), y(r2)-y(r1), Lx, Ly);

        nbond = nbond + 1;
        BondsRows(nbond,:) = [r1, r2, L];

        deg(r1) = deg(r1) + 1;
        deg(r2) = deg(r2) + 1;
        adj(r1,r2) = 1;
        adj(r2,r1) = 1;

        no_progress = 0;
    end
    obj.log.print('   Mono/random: placed %d bonds in %4.4f sec\n', nbond, toc);

    BondsRows = BondsRows(1:nbond,:);

    skip_prune = ~isempty(PreBonds);
    [AtomsOut, BondsOut] = finalize_network(obj, Atoms, BondsRows, ids, x, y, ...
        Max_peratom_bond, min_keep, isPeriodic, Lx, Ly, skip_prune);

end


% =========================================================================
% LATTICE MONO
% =========================================================================
function [Atoms, Bonds] = connect_lattice_mono(obj, Atoms, LatticeData)

    if isempty(LatticeData) || ~isstruct(LatticeData) || ~isfield(LatticeData,'idx_map')
        error('AddBondsMono: lattice geometry requires LatticeData.idx_map.');
    end

    idx_map = LatticeData.idx_map;
    Ny = size(idx_map,1);
    Nx = size(idx_map,2);

    Natoms = size(Atoms,1);
    maxBonds = 3 * Natoms;

    i_list  = zeros(maxBonds,1);
    j_list  = zeros(maxBonds,1);
    L0_list = zeros(maxBonds,1);
    nb = 0;

    for iy = 1:Ny
        for ix = 1:Nx

            i_atom = idx_map(iy, ix);
            if i_atom == 0
                continue;
            end

            % Horizontal forward neighbor
            if ix < Nx
                j_atom = idx_map(iy, ix+1);
                if j_atom > 0
                    nb = nb + 1;
                    i_list(nb) = i_atom;
                    j_list(nb) = j_atom;

                    dx = Atoms(j_atom,3) - Atoms(i_atom,3);   % X at col 3
                    dy = Atoms(j_atom,4) - Atoms(i_atom,4);   % Y at col 4
                    L0_list(nb) = sqrt(dx*dx + dy*dy);
                end
            end

            % Two forward diagonal neighbors
            if iy < Ny
                if mod(iy,2) == 1
                    j1 = idx_map(iy+1, ix);
                    j2 = 0;
                    if ix < Nx
                        j2 = idx_map(iy+1, ix+1);
                    end
                else
                    j1 = idx_map(iy+1, ix);
                    j2 = 0;
                    if ix > 1
                        j2 = idx_map(iy+1, ix-1);
                    end
                end

                if j1 > 0
                    nb = nb + 1;
                    i_list(nb) = i_atom;
                    j_list(nb) = j1;

                    dx = Atoms(j1,3) - Atoms(i_atom,3);
                    dy = Atoms(j1,4) - Atoms(i_atom,4);
                    L0_list(nb) = sqrt(dx*dx + dy*dy);
                end

                if j2 > 0
                    nb = nb + 1;
                    i_list(nb) = i_atom;
                    j_list(nb) = j2;

                    dx = Atoms(j2,3) - Atoms(i_atom,3);
                    dy = Atoms(j2,4) - Atoms(i_atom,4);
                    L0_list(nb) = sqrt(dx*dx + dy*dy);
                end
            end
        end
    end

    i_list  = i_list(1:nb);
    j_list  = j_list(1:nb);
    L0_list = L0_list(1:nb);

    Bonds = zeros(nb,5);
    Bonds(:,1) = (1:nb).';
    Bonds(:,2) = i_list;
    Bonds(:,3) = j_list;
    Bonds(:,4) = L0_list;
    Bonds(:,5) = 0;   % type left for AssignMultiType

    % Rebuild degree / neighbor list in the new column layout
    Max_peratom_bond = obj.peratom.Max_peratom_bond;
    needed_cols = 6 + Max_peratom_bond;
    if size(Atoms,2) < needed_cols
        Atoms(:, size(Atoms,2)+1 : needed_cols) = 0;
    end

    Atoms(:,6) = 0;                            % degree at col 6
    Atoms(:,7:6+Max_peratom_bond) = 0;         % neighbor slots at 7..

    for k = 1:size(Bonds,1)
        ii = Bonds(k,2);
        jj = Bonds(k,3);

        Atoms(ii,6) = Atoms(ii,6) + 1;
        if Atoms(ii,6) <= Max_peratom_bond
            Atoms(ii, 6 + Atoms(ii,6)) = jj;
        end

        Atoms(jj,6) = Atoms(jj,6) + 1;
        if Atoms(jj,6) <= Max_peratom_bond
            Atoms(jj, 6 + Atoms(jj,6)) = ii;
        end
    end

    obj.log.print('   Mono/lattice: placed %d bonds\n', size(Bonds,1));

end


% =========================================================================
% HELPERS
% =========================================================================
function neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic)

    Cx = cx(r1);
    Cy = cy(r1);
    neigh = [];

    if isPeriodic
        for dxCell = -1:1
            ix = Cx + dxCell;
            if ix < 1
                ix = nx;
            elseif ix > nx
                ix = 1;
            end

            for dyCell = -1:1
                iy = Cy + dyCell;
                if iy < 1
                    iy = ny;
                elseif iy > ny
                    iy = 1;
                end

                if ~isempty(Cells{ix,iy})
                    neigh = [neigh, Cells{ix,iy}]; %#ok<AGROW>
                end
            end
        end
    else
        for dxCell = -1:1
            ix = Cx + dxCell;
            if ix < 1 || ix > nx
                continue;
            end

            for dyCell = -1:1
                iy = Cy + dyCell;
                if iy < 1 || iy > ny
                    continue;
                end

                if ~isempty(Cells{ix,iy})
                    neigh = [neigh, Cells{ix,iy}]; %#ok<AGROW>
                end
            end
        end
    end

end

function d = minimum_image(isPeriodic, dx, dy, Lx, Ly)

    if ~isPeriodic
        d = sqrt(dx.^2 + dy.^2);
        return;
    end

    dxp = dx - Lx .* round(dx ./ Lx);
    dyp = dy - Ly .* round(dy ./ Ly);
    d = sqrt(dxp.^2 + dyp.^2);

end

function [AtomsOut, BondsOut] = finalize_network(obj, Atoms, BondsRows, ids, x, y, ...
    Max_peratom_bond, min_keep, isPeriodic, Lx, Ly, skip_prune)
% Finalizes bonds for the random-mono path.
%   - If skip_prune == true (multi-phase build, e.g. bottle-brush crosslinks),
%     internal low-degree pruning is bypassed. The AddBonds dispatcher
%     rebuilds neighbor lists on the combined rod+crosslink bond list, and
%     CleanupNetwork does the global pruning later.
%   - Bond type (col 5) is left 0; AssignMultiType assigns types.

    natom = size(Atoms,1);
    pruned_atom_mask = false(natom,1);

    if ~isempty(BondsRows) && ~skip_prune
        changed = true;
        while changed
            deg = zeros(natom,1);
            for k = 1:size(BondsRows,1)
                deg(BondsRows(k,1)) = deg(BondsRows(k,1)) + 1;
                deg(BondsRows(k,2)) = deg(BondsRows(k,2)) + 1;
            end

            to_del = find(deg <= min_keep);
            if isempty(to_del)
                changed = false;
                break;
            end

            kill = false(size(BondsRows,1),1);
            mark = false(natom,1);
            mark(to_del) = true;

            for k = 1:size(BondsRows,1)
                if mark(BondsRows(k,1)) || mark(BondsRows(k,2))
                    kill(k) = true;
                end
            end

            if any(kill)
                BondsRows = BondsRows(~kill,:);
                pruned_atom_mask(to_del) = true;
                changed = true;
            else
                changed = false;
            end
        end
    end

    if any(pruned_atom_mask)
        Atoms = Atoms(~pruned_atom_mask, :);
        ids   = Atoms(:,1);
        x     = Atoms(:,3);
        y     = Atoms(:,4);
        natom = size(Atoms,1);

        old2new_row = zeros(numel(pruned_atom_mask),1,'int32');
        old2new_row(~pruned_atom_mask) = int32(1:natom);

        if ~isempty(BondsRows)
            keepBond = ~pruned_atom_mask(BondsRows(:,1)) & ~pruned_atom_mask(BondsRows(:,2));
            BondsRows = BondsRows(keepBond,:);
            BondsRows(:,1) = old2new_row(BondsRows(:,1));
            BondsRows(:,2) = old2new_row(BondsRows(:,2));
        end
    end

    if isempty(Atoms) || isempty(BondsRows)
        AtomsOut = Atoms;
        BondsOut = zeros(0,5);
        if ~skip_prune
            obj.log.print('   Pruned all atoms/bonds\n');
        end
        return;
    end

    for k = 1:size(BondsRows,1)
        r1 = BondsRows(k,1);
        r2 = BondsRows(k,2);
        BondsRows(k,3) = minimum_image(isPeriodic, x(r2)-x(r1), y(r2)-y(r1), Lx, Ly);
    end

    nb = size(BondsRows,1);
    BondsOut = zeros(nb,5);
    for k = 1:nb
        r1 = BondsRows(k,1);
        r2 = BondsRows(k,2);
        BondsOut(k,:) = [k, ids(r1), ids(r2), BondsRows(k,3), 0];  % type left 0
    end

    oldIDs = ids;
    newIDs = (1:natom).';
    Atoms(:,1) = newIDs;

    [tfI,locI] = ismember(BondsOut(:,2), oldIDs);
    [tfJ,locJ] = ismember(BondsOut(:,3), oldIDs);
    BondsOut(tfI,2) = newIDs(locI(tfI));
    BondsOut(tfJ,3) = newIDs(locJ(tfJ));

    if ~isempty(BondsOut)
        BondsOut(:,1) = (1:size(BondsOut,1)).';
    end

    % Rebuild per-atom neighbor lists in the new column layout
    needed_cols = 6 + Max_peratom_bond;
    if size(Atoms,2) < needed_cols
        Atoms(:, size(Atoms,2)+1 : needed_cols) = 0;
    end

    Atoms(:,6) = 0;
    Atoms(:,7:6+Max_peratom_bond) = 0;

    for k = 1:size(BondsOut,1)
        ii = BondsOut(k,2);
        jj = BondsOut(k,3);

        nb1 = Atoms(ii,6) + 1;
        if nb1 <= Max_peratom_bond
            Atoms(ii,6) = nb1;
            Atoms(ii, 6 + nb1) = jj;
        end

        nb2 = Atoms(jj,6) + 1;
        if nb2 <= Max_peratom_bond
            Atoms(jj,6) = nb2;
            Atoms(jj, 6 + nb2) = ii;
        end
    end

    AtomsOut = Atoms;

end
