function [Atoms, Bonds] = AddBondsBimodal(obj, Atoms, LatticeData, PreBonds)
% -------------------------------------------------------------------------
% AddBondsBimodal
%
% Bimodal bond-topology generation:
%   - random / hex_lattice / bottle_brush all route to the same generalized
%     bimodal distance-window connector.
%
% Atoms column layout:
%   [ID | molID | X | Y | Z | deg | nbr_1 ... nbr_{Max_peratom_bond}]
%
% INPUT:
%   obj         : network object
%   Atoms       : atom array
%   LatticeData : unused (kept for dispatcher symmetry)
%   PreBonds    : optional pre-existing bonds; seeds adjacency so the long
%                 and short bond passes don't duplicate rod bonds. When
%                 non-empty, internal pruning is skipped.
%
% OUTPUT:
%   Atoms : updated atom array
%   Bonds : [bondID | atomID_i | atomID_j | L0 | type=0]
%           NOTE: bond type is LEFT ZERO. The connector uses internal
%           long/short labels only for its own bookkeeping and clears the
%           column before returning - AssignMultiType decides real types.
% -------------------------------------------------------------------------

    if nargin < 4
        PreBonds = [];
    end

    [Atoms, Bonds] = connect_bimodal_general(obj, Atoms, PreBonds);

    %#ok<NASGU>

end


% =========================================================================
% GENERAL BIMODAL CONNECTOR
% =========================================================================
function [Atoms, Bonds] = connect_bimodal_general(obj, Atoms, PreBonds)

    natom            = size(Atoms,1);
    Max_bond         = obj.domain.Max_bond;
    Max_peratom_bond = obj.peratom.Max_peratom_bond;
    xlo = obj.domain.xlo; xhi = obj.domain.xhi;
    ylo = obj.domain.ylo; yhi = obj.domain.yhi;

    Lx = xhi - xlo;
    Ly = yhi - ylo;
    domain_diag = hypot(Lx, Ly);

    global_limit = obj.domain.bond_global_try_limit;
    stall_limit  = obj.domain.max_attempts_without_progress;
    min_keep     = obj.peratom.min_degree_keep;
    dmin         = obj.domain.min_node_sep;
    epsr         = 1e-9;

    isPeriodic = strcmpi(obj.domain.boundary, 'periodic');
    skip_prune = ~isempty(PreBonds);

    b = obj.domain.b;

    bi = obj.architecture.strand_typology.bimodal;

    N1   = bi.mean_1;
    N2   = bi.mean_2;
    sig1 = bi.std_1;
    sig2 = bi.std_2;
    sigr1 = bi.stdR_1;
    sigr2 = bi.stdR_2;
    lam1 = bi.lam_1;
    lam2 = bi.lam_2;

    useProb    = strcmpi(bi.height_mode, 'prob');
    useManual  = strcmpi(bi.bin_window_method, 'manual');
    useMixed   = strcmpi(bi.manual_dev_type, 'mixed') || ...
                 strcmpi(bi.manual_dev_type, 'both');
    long_first = bi.long_first;

    double_network = bi.double_network_flag;
    autoN1 = bi.auto_1_flag;
    autoN2 = bi.auto_2_flag;

    if double_network
        alpha = bi.alpha;
        if isempty(alpha)
            alpha = 2.5;
        end
        f_sparse = min(1, 1/(alpha^2));
    else
        alpha = 1.0;
        f_sparse = 1.0;
    end

    if useProb
        P2 = bi.height_prob;
        target_N2 = max(0, min(Max_bond, round(P2 * Max_bond)));
    else
        target_N2 = bi.height_count;
        target_N2 = max(0, min(Max_bond, target_N2));
    end

    z1_min = 3;
    targetC1 = 8;
    targetC2 = 12;
    Cmin = 4;
    Cmax = 24;

    if lam1 < 0 || lam1 > 1
        lam1 = 1 / sqrt(max(N1,1));
    end
    if lam2 < 0 || lam2 > 1
        lam2 = 1 / sqrt(max(N2,1));
    end

    r1_avg = lam1 * b * N1;
    r_min_allowed = max(dmin * 1.2, b * 0.5);
    if r1_avg < r_min_allowed
        r1_avg = r_min_allowed;
    end

    if double_network
        r2_avg = alpha * r1_avg;
    else
        r2_avg = lam2 * b * N2;
    end

    if double_network && r2_avg < r_min_allowed
        r2_avg = r1_avg;
    end

    if (r2_avg < 1.8 * r1_avg) && (~double_network)
        r2_avg = 1.8 * r1_avg;
    end

    r2_max_allowed = 0.4 * domain_diag;
    if r2_avg > r2_max_allowed
        r2_avg = r2_max_allowed;
    end

    A = (xhi-xlo) * (yhi-ylo);
    rho_all = natom / max(A, epsr);

    if useManual
        if useMixed
            dr1 = 2.355 * sigr1;
            dr2 = 2.355 * sigr2;
        else
            dr1 = lam1 * b * (2.355 * sig1);
            dr2 = lam2 * b * (2.355 * sig2);
        end
    else
        dr1 = targetC1 / max(2*pi*rho_all*max(r1_avg,epsr), epsr);

        Nsparse_est = max(1, round(f_sparse * natom));
        if double_network
            rho_long = Nsparse_est / max(A, epsr);
        else
            rho_long = rho_all;
        end
        dr2 = targetC2 / max(2*pi*rho_long*max(r2_avg,epsr), epsr);
    end

    r1_lower = max(r1_avg - dr1, dmin + epsr);
    r1_upper = r1_avg + dr1;

    gap = 0.1 * r1_avg;
    r2_lower = r2_avg - dr2;
    r2_upper = r2_avg + dr2;

    if autoN1
        N1old = N1;
        R1AVG = 0.5*(r1_upper-r1_lower) + r1_lower;
        N1 = R1AVG/(lam1*b);
        obj.architecture.strand_typology.bimodal.mean_1 = N1;
        obj.log.print('   Auto N1: adjusted N1 from %.0f to %.0f\n', N1old, N1);
    end

    if autoN2
        N2old = N2;
        R2AVG = 0.5*(r2_upper-r2_lower) + r2_lower;
        N2 = 1.4*R2AVG/(lam2*b);
        obj.architecture.strand_typology.bimodal.mean_2 = N2;
        obj.log.print('   Auto N2: adjusted N2 from %.0f to %.0f\n', N2old, N2);
    end

    ids = Atoms(:,1);
    mol_ids = Atoms(:,2);
    sameRod = strcmpi(obj.architecture.geometry, 'bottle_brush');
    x   = Atoms(:,3);   % X at col 3
    y   = Atoms(:,4);   % Y at col 4

    % Live per-node adjacency set used by exclude_existing_any.
    adj_set = cell(natom, 1);

    % Seed adjacency from PreBonds so duplicate bonds aren't generated.
    % Atom IDs in PreBonds equal row indices at this point.
    if ~isempty(PreBonds)
        for k = 1:size(PreBonds, 1)
            r1 = PreBonds(k, 2);
            r2 = PreBonds(k, 3);
            adj_set{r1}(end+1) = r2; %#ok<AGROW>
            adj_set{r2}(end+1) = r1; %#ok<AGROW>
        end
    end

    if double_network
        avg_nn_spacing = sqrt((xhi-xlo)*(yhi-ylo)/natom);
        target_spacing = alpha * avg_nn_spacing;
        if alpha == 1
            target_spacing = 10;
        end
        isSparse = pick_uniform_sparse_nodes(x, y, f_sparse, target_spacing);
        sparse_idx = find(isSparse);
    else
        isSparse = true(natom,1);
        sparse_idx = 1:natom;
    end

    rmax = max(r1_upper, r2_upper);
    if rmax <= 0
        rmax = max(xhi-xlo, yhi-ylo);
    end
    cellSize = rmax;

    nx = max(1, ceil((xhi-xlo)/cellSize));
    ny = max(1, ceil((yhi-ylo)/cellSize));

    cx = floor((x - xlo)/cellSize) + 1;
    cy = floor((y - ylo)/cellSize) + 1;
    cx = max(1, min(nx, cx));
    cy = max(1, min(ny, cy));

    Cells = cell(nx, ny);
    for i = 1:natom
        Cells{cx(i), cy(i)}(end+1) = i; %#ok<AGROW>
    end

    deg1 = zeros(natom,1);
    deg2 = zeros(natom,1);

    Btmp = zeros(Max_bond, 5);   % [bid r1 r2 L internal_label(1/2)]
    nbond = 0;
    countType2 = 0;

    g_mul1 = 1.0;
    g_mul2 = 1.0;

    if long_first
        % ---- long bonds first ----
        ntries = 0;
        no_progress = 0;

        while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
            if countType2 >= target_N2
                break;
            end

            ntries = ntries + 1;

            if double_network
                r1 = sparse_idx(randi(numel(sparse_idx)));
            else
                r1 = randi(natom);
            end

            if deg2(r1) >= Max_peratom_bond
                no_progress = no_progress + 1;
                continue;
            end

            if double_network
                r2lo = r2_lower;
                r2hi = r2_upper;
            else
                dr2_pick = dr2 * g_mul2;
                r2lo = max(r2_lower - (dr2-dr2_pick), r1_upper + gap);
                r2hi = r2_upper + (dr2_pick-dr2);
            end

            neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            if double_network
                neigh = neigh(isSparse(neigh));
                if isempty(neigh)
                    no_progress = no_progress + 1;
                    continue;
                end
            end

            neigh(neigh==r1) = [];
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            neigh = neigh(deg2(neigh) < Max_peratom_bond);
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            dxv = x(neigh) - x(r1);
            dyv = y(neigh) - y(r1);
            d = minimum_image(isPeriodic, dxv, dyv, Lx, Ly);

            in2 = (d >= r2lo) & (d <= r2hi);
            cand2 = neigh(in2);
            cand2 = exclude_existing_any(cand2, r1, adj_set);
            if sameRod
                cand2 = cand2(mol_ids(cand2) ~= mol_ids(r1));
            end

            if ~useManual
                C = numel(cand2);
                if C < Cmin
                    g_mul2 = min(2.0, g_mul2*1.15);
                    no_progress = no_progress + 1;
                    continue;
                elseif C > Cmax
                    g_mul2 = max(0.5, g_mul2*0.85);
                else
                    g_mul2 = 1.0;
                end
            end

            if isempty(cand2)
                no_progress = no_progress + 1;
                continue;
            end

            [~,ord] = sort(deg2(cand2), 'ascend');
            cand2 = cand2(ord);
            r2 = cand2(min(numel(cand2), randi(min(5,numel(cand2)))));

            L = minimum_image(isPeriodic, x(r2)-x(r1), y(r2)-y(r1), Lx, Ly);

            nbond = nbond + 1;
            Btmp(nbond,:) = [nbond, r1, r2, L, 2];
            countType2 = countType2 + 1;

            deg2(r1) = deg2(r1) + 1;
            deg2(r2) = deg2(r2) + 1;

            adj_set{r1}(end+1) = r2; %#ok<AGROW>
            adj_set{r2}(end+1) = r1; %#ok<AGROW>

            no_progress = 0;
        end

        % ---- short bonds after ----
        ntries = 0;
        no_progress = 0;

        while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)

            ntries = ntries + 1;

            r1 = randi(natom);
            if deg1(r1) >= Max_peratom_bond
                no_progress = no_progress + 1;
                continue;
            end

            dr1_pick = dr1 * g_mul1;
            rlo = max(r1_lower - (dr1-dr1_pick), dmin + epsr);
            rhi = r1_upper + (dr1_pick-dr1);

            neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            neigh(neigh==r1) = [];
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            neigh = neigh(deg1(neigh) < Max_peratom_bond);
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            dxv = x(neigh) - x(r1);
            dyv = y(neigh) - y(r1);
            d = minimum_image(isPeriodic, dxv, dyv, Lx, Ly);

            in1 = (d >= rlo) & (d <= rhi);
            cand = neigh(in1);
            cand = exclude_existing_any(cand, r1, adj_set);
            if sameRod
                cand = cand(mol_ids(cand) ~= mol_ids(r1));
            end

            if ~useManual
                C = numel(cand);
                if C < Cmin
                    g_mul1 = min(2.0, g_mul1*1.15);
                    no_progress = no_progress + 1;
                    continue;
                elseif C > Cmax
                    g_mul1 = max(0.5, g_mul1*0.85);
                else
                    g_mul1 = 1.0;
                end
            end

            if isempty(cand)
                no_progress = no_progress + 1;
                continue;
            end

            r2 = cand(randi(numel(cand)));
            L = minimum_image(isPeriodic, x(r2)-x(r1), y(r2)-y(r1), Lx, Ly);

            nbond = nbond + 1;
            Btmp(nbond,:) = [nbond, r1, r2, L, 1];

            deg1(r1) = deg1(r1) + 1;
            deg1(r2) = deg1(r2) + 1;

            adj_set{r1}(end+1) = r2; %#ok<AGROW>
            adj_set{r2}(end+1) = r1; %#ok<AGROW>

            no_progress = 0;
        end

    else
        % ---- short scaffold first ----
        ntries = 0;
        no_progress = 0;
        degTot = zeros(natom,1);

        while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
            if mean(degTot >= z1_min) > 0.95
                break;
            end

            ntries = ntries + 1;

            r1 = randi(natom);
            if degTot(r1) >= Max_peratom_bond
                no_progress = no_progress + 1;
                continue;
            end

            dr1_pick = dr1 * g_mul1;
            rlo = max(r1_lower - (dr1-dr1_pick), dmin + epsr);
            rhi = r1_upper + (dr1_pick-dr1);

            neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            neigh(neigh==r1) = [];
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            neigh = neigh(degTot(neigh) < Max_peratom_bond);
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            dxv = x(neigh) - x(r1);
            dyv = y(neigh) - y(r1);
            d = minimum_image(isPeriodic, dxv, dyv, Lx, Ly);

            in1 = (d >= rlo) & (d <= rhi);
            cand = neigh(in1);
            cand = exclude_existing_any(cand, r1, adj_set);
            if sameRod
                cand = cand(mol_ids(cand) ~= mol_ids(r1));
            end

            if ~useManual
                C = numel(cand);
                if C < Cmin
                    g_mul1 = min(2.0, g_mul1*1.15);
                    no_progress = no_progress + 1;
                    continue;
                elseif C > Cmax
                    g_mul1 = max(0.5, g_mul1*0.85);
                else
                    g_mul1 = 1.0;
                end
            end

            if isempty(cand)
                no_progress = no_progress + 1;
                continue;
            end

            r2 = cand(randi(numel(cand)));
            L = minimum_image(isPeriodic, x(r2)-x(r1), y(r2)-y(r1), Lx, Ly);

            nbond = nbond + 1;
            Btmp(nbond,:) = [nbond, r1, r2, L, 1];

            deg1(r1) = deg1(r1) + 1;
            deg1(r2) = deg1(r2) + 1;
            degTot(r1) = degTot(r1) + 1;
            degTot(r2) = degTot(r2) + 1;

            adj_set{r1}(end+1) = r2; %#ok<AGROW>
            adj_set{r2}(end+1) = r1; %#ok<AGROW>

            no_progress = 0;
        end

        % ---- long bonds second ----
        ntries = 0;
        no_progress = 0;

        while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
            if (~useProb) && (countType2 >= target_N2)
                break;
            end
            if useProb && (countType2 > target_N2*1.1)
                break;
            end

            ntries = ntries + 1;

            if double_network
                r1 = sparse_idx(randi(numel(sparse_idx)));
            else
                r1 = randi(natom);
            end

            if degTot(r1) >= Max_peratom_bond
                no_progress = no_progress + 1;
                continue;
            end

            if double_network
                r2lo = r2_lower;
                r2hi = r2_upper;
            else
                dr2_pick = dr2 * g_mul2;
                r2lo = max(r2_lower - (dr2-dr2_pick), r1_upper + gap);
                r2hi = r2_upper + (dr2_pick-dr2);
            end

            neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            if double_network
                neigh = neigh(isSparse(neigh));
                if isempty(neigh)
                    no_progress = no_progress + 1;
                    continue;
                end
            end

            neigh(neigh==r1) = [];
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            neigh = neigh(degTot(neigh) < Max_peratom_bond);
            if isempty(neigh)
                no_progress = no_progress + 1;
                continue;
            end

            dxv = x(neigh) - x(r1);
            dyv = y(neigh) - y(r1);
            d = minimum_image(isPeriodic, dxv, dyv, Lx, Ly);

            in2 = (d >= r2lo) & (d <= r2hi);
            cand2 = neigh(in2);
            cand2 = exclude_existing_any(cand2, r1, adj_set);
            if sameRod
                cand2 = cand2(mol_ids(cand2) ~= mol_ids(r1));
            end

            if ~useManual
                C = numel(cand2);
                if C < Cmin
                    g_mul2 = min(2.0, g_mul2*1.15);
                    no_progress = no_progress + 1;
                    continue;
                elseif C > Cmax
                    g_mul2 = max(0.5, g_mul2*0.85);
                else
                    g_mul2 = 1.0;
                end
            end

            if isempty(cand2)
                no_progress = no_progress + 1;
                continue;
            end

            [~,ord] = sort(degTot(cand2), 'ascend');
            cand2 = cand2(ord);
            r2 = cand2(min(numel(cand2), randi(min(5,numel(cand2)))));

            L = minimum_image(isPeriodic, x(r2)-x(r1), y(r2)-y(r1), Lx, Ly);

            nbond = nbond + 1;
            Btmp(nbond,:) = [nbond, r1, r2, L, 2];
            countType2 = countType2 + 1;

            deg2(r1) = deg2(r1) + 1;
            deg2(r2) = deg2(r2) + 1;
            degTot(r1) = degTot(r1) + 1;
            degTot(r2) = degTot(r2) + 1;

            adj_set{r1}(end+1) = r2; %#ok<AGROW>
            adj_set{r2}(end+1) = r1; %#ok<AGROW>

            no_progress = 0;
        end
    end

    Btmp = Btmp(1:nbond,:);

    [Atoms, Bonds] = finalize_bimodal_network(Atoms, Btmp, ids, Max_peratom_bond, ...
        min_keep, isPeriodic, Lx, Ly, skip_prune);

    % Bond type column (5) is set to 0: AssignMultiType is authoritative.
    % The connector's internal long/short label (in Btmp col 5) is discarded
    % by the finalize helper, which writes col 5 = 0.
    obj.log.print('   Bimodal/%s: placed %d bonds (types left to AssignMultiType)\n', ...
        lower(obj.architecture.geometry), size(Bonds,1));

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

function cand = exclude_existing_any(cand, r1, adj_set)

    if isempty(cand) || isempty(adj_set{r1})
        return;
    end

    already = adj_set{r1};
    cand    = cand(~ismember(cand, already));

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

function isSparse = pick_uniform_sparse_nodes(x, y, f_sparse, r_spacing)

    natom          = numel(x);
    Nsparse_target = round(f_sparse * natom);
    isSparse       = false(natom, 1);

    idx_all = randperm(natom);
    picked  = [];

    for k = 1:natom
        i = idx_all(k);

        if isempty(picked)
            picked(end+1) = i; %#ok<AGROW>
            continue;
        end

        dx = x(picked) - x(i);
        dy = y(picked) - y(i);
        if all(sqrt(dx.^2 + dy.^2) >= r_spacing)
            picked(end+1) = i; %#ok<AGROW>
        end

        if numel(picked) >= Nsparse_target
            break;
        end
    end

    isSparse(picked) = true;

    if numel(picked) < Nsparse_target
        fprintf(['   [pick_uniform_sparse_nodes] Warning: placed %d / %d sparse nodes.\n' ...
                 '   Domain may be too small or r_spacing too large for requested f_sparse.\n'], ...
                numel(picked), Nsparse_target);
    end

end

function [Atoms, Bonds] = finalize_bimodal_network(Atoms, Btmp, ids, ...
    Max_peratom_bond, min_keep, isPeriodic, Lx, Ly, skip_prune)
% New-layout finalize helper. Bond type col (5) is written as 0 regardless
% of the internal label carried in Btmp: AssignMultiType sets the real type
% later.

    natom = size(Atoms,1);

    if ~isempty(Btmp) && (min_keep > 0) && ~skip_prune
        changed = true;
        while changed
            deg_tot = zeros(natom,1);
            for k = 1:size(Btmp,1)
                deg_tot(Btmp(k,2)) = deg_tot(Btmp(k,2)) + 1;
                deg_tot(Btmp(k,3)) = deg_tot(Btmp(k,3)) + 1;
            end

            to_del = find(deg_tot <= min_keep);
            if isempty(to_del)
                changed = false;
                break;
            end

            kill = false(size(Btmp,1),1);
            mark = false(natom,1);
            mark(to_del) = true;

            for k = 1:size(Btmp,1)
                if mark(Btmp(k,2)) || mark(Btmp(k,3))
                    kill(k) = true;
                end
            end

            if any(kill)
                Btmp = Btmp(~kill,:);
                changed = true;
            else
                changed = false;
            end
        end
    end

    Bonds = zeros(size(Btmp,1), 5);
    for k = 1:size(Btmp,1)
        r1 = Btmp(k,2);
        r2 = Btmp(k,3);
        Bonds(k,:) = [k, ids(r1), ids(r2), Btmp(k,4), 0];
    end

    % Second pruning pass operates in atom-ID space, iteratively removing
    % low-degree atoms until converged. Bypassed when skip_prune is true.
    if (min_keep > 0) && ~skip_prune
        changed = true;
        while changed
            changed = false;

            if isempty(Bonds)
                Atoms = [];
                break;
            end

            id2row = containers.Map('KeyType','int64','ValueType','int32');
            for r = 1:size(Atoms,1)
                id2row(int64(Atoms(r,1))) = int32(r);
            end

            nA = size(Atoms,1);
            deg = zeros(nA,1);
            ri = zeros(size(Bonds,1),1,'int32');
            rj = ri;

            for k = 1:size(Bonds,1)
                ri(k) = id2row(int64(Bonds(k,2)));
                rj(k) = id2row(int64(Bonds(k,3)));
            end

            idx = [ri; rj];
            deg = deg + accumarray(double(idx), ones(numel(idx),1), [nA 1], @sum, 0);

            delMask = (deg <= min_keep);
            if ~any(delMask)
                break;
            end

            delIDs = Atoms(delMask,1);

            keepBond = ~ismember(Bonds(:,2), delIDs) & ~ismember(Bonds(:,3), delIDs);
            if any(~keepBond)
                Bonds = Bonds(keepBond,:);
                changed = true;
            end

            keepAtom = ~delMask;
            if any(~keepAtom)
                Atoms = Atoms(keepAtom,:);
                changed = true;
            end
        end
    end

    if isempty(Atoms) || isempty(Bonds)
        Atoms = zeros(0, 6 + Max_peratom_bond);
        Bonds = zeros(0, 5);
        return;
    end

    oldIDs = Atoms(:,1);
    newIDs = (1:size(Atoms,1)).';
    Atoms(:,1) = newIDs;

    [tfI,locI] = ismember(Bonds(:,2), oldIDs);
    [tfJ,locJ] = ismember(Bonds(:,3), oldIDs);
    Bonds(tfI,2) = newIDs(locI(tfI));
    Bonds(tfJ,3) = newIDs(locJ(tfJ));
    Bonds(:,1) = (1:size(Bonds,1)).';

    % Recompute bond lengths under current (possibly renumbered) atom positions
    for k = 1:size(Bonds,1)
        i = Bonds(k,2);
        j = Bonds(k,3);
        Bonds(k,4) = minimum_image(isPeriodic, ...
            Atoms(j,3)-Atoms(i,3), Atoms(j,4)-Atoms(i,4), Lx, Ly);
    end

    % Rebuild neighbor lists in the new column layout
    needed_cols = 6 + Max_peratom_bond;
    if size(Atoms,2) < needed_cols
        Atoms(:, size(Atoms,2)+1 : needed_cols) = 0;
    end

    Atoms(:,6) = 0;
    Atoms(:,7:6+Max_peratom_bond) = 0;

    for k = 1:size(Bonds,1)
        ii = Bonds(k,2);
        jj = Bonds(k,3);

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

end
