function [Atoms,Bonds,options] = NetworkGenConnectNodesBimodal(Domain,Atoms,options,advancedOptions)
% NetworkGenConnectNodesBimodal (adaptive, grid-accelerated)
% Adds 'long_first' mode: place type-2 bonds first (fixed N2_number),
% then place type-1 bonds as if type-2 bonds don't consume capacity.
% Adds 'double_network' mode (in the same spirit as long_first):
%   - pick a sparse subset of nodes (2–3x larger spacing),
%   - place type-2 bonds only among that sparse subset,
%   - (optional) size long-bond window using the sparse density.
%
% Returns Bonds(:,5) = 1 or 2 for type.

% ---------- Unpack domain & defaults ----------
natom            = size(Atoms,1);
Max_bond         = Domain.Max_bond;
Max_peratom_bond = Domain.Max_peratom_bond;
xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;

Lx = xhi - xlo; Ly = yhi - ylo;
domain_diag = hypot(Lx, Ly);

if isfield(Domain,'bond_global_try_limit'), global_limit = Domain.bond_global_try_limit; else, global_limit = 5e7; end
if isfield(Domain,'max_attempts_without_progress'), stall_limit = Domain.max_attempts_without_progress; else, stall_limit = 5e5; end
if isfield(Domain,'min_degree_keep'), min_keep = Domain.min_degree_keep; else, min_keep = 1; end
if isfield(Domain,'min_node_sep'), dmin = Domain.min_node_sep; else, dmin = 0; end
epsr = 1e-9;

% ---------- Options ----------
b    = options.b;
N1   = options.bimodal.N1;
N2   = options.bimodal.N2;
sig1 = options.bimodal.std1;
sig2 = options.bimodal.std2;
sigr1 = options.bimodal.stdr1;
sigr2 = options.bimodal.stdr2;
lam1 = options.bimodal.lam1;
lam2 = options.bimodal.lam2;

% Read mode options
useProb     = strcmpi(options.bimodal.distribution_height_mode,'prob');
useManual   = strcmpi(options.bimodal.bin_window_method,'manual');
useMixed    = strcmpi(options.bimodal.deviation_type,'mixed');
long_first  = isfield(options.bimodal,'long_first') && options.bimodal.long_first;
isPeriodic  =  strcmpi(options.boundary_box,'Periodic');

% --- Double-network flag & params (minimal style like long_first) ---
double_network = isfield(options,'double_network') && (options.double_network.flag);
autoN2 = isfield(options,'double_network') && (options.double_network.autoN2);

if double_network
    dn = options.double_network;
    % Either specify alpha (spacing multiplier wrt small mesh) or sparse_fraction directly
    if isfield(dn,'sparse_fraction') && ~isempty(dn.sparse_fraction)
        f_sparse = max(0,min(1,dn.sparse_fraction));
    else
        if ~isfield(dn,'alpha') || isempty(dn.alpha), dn.alpha = 2.5; end % default spacing ~2.5x
        f_sparse = min(1, 1/(dn.alpha^2)); % 2D: spacing ~ 1/sqrt(rho)
    end
else
    f_sparse = 1.0; % degenerate: long pass can touch all nodes
end

if useProb
    P2 = options.bimodal.P;
    target_N2 = max(0, min(Max_bond, round(P2*Max_bond)));
else
    target_N2 = options.bimodal.N2_number;
    target_N2 = max(0, min(Max_bond, target_N2));
end

if isfield(options.bimodal,'z1_min'), z1_min = options.bimodal.z1_min; else, z1_min = 3; end
if isfield(options.bimodal,'targetC1'), targetC1 = options.bimodal.targetC1; else, targetC1 = 8; end
if isfield(options.bimodal,'targetC2'), targetC2 = options.bimodal.targetC2; else, targetC2 = 12; end
Cmin = 4; Cmax = 24;

% ---------- Apply checks to the prestretch variables ----------
if lam1 <0 || lam1 >1
    warning('lam1=%.3g out of range [0,1]; reverting to default: 1/sqrt(N1)=%.3g.', lam1, 1/sqrt(N1));
    lam1 = 1/sqrt(N1);
end
if lam2 <0 || lam2 >1
    
    warning('lam2=%.3g out of range [0,1]; reverting to default: 1/sqrt(N2)=%.3g.', lam2, 1/sqrt(N2));
    lam2 = 1/sqrt(N2);
end

% ---------- Refine r1_avg and r2_avg based on geometry ----------
r1_avg = lam1*b*N1;

%check r1_avg vs min allowed spacing
r_min_allowed = max(dmin * 1.2, b * 0.5);
if r1_avg < r_min_allowed
    warning('r1_avg=%.3g below min separation %.3g; lifting to %.3g.', r1_avg, dmin, r_min_allowed);
    r1_avg = r_min_allowed;
end

% if double network then scale
if double_network
    r2_avg = dn.alpha*r1_avg;
else
    r2_avg = lam2*b*N2;
end

if double_network && r2_avg < r_min_allowed
    r2_avg = r1_avg;
end

if r2_avg < 1.8 * r1_avg && ~double_network
    warning('r2_avg too close to r1_avg (%.3g vs %.3g). Adjusting r2_avg.', r2_avg, r1_avg);
    r2_avg = 1.8 * r1_avg;
end

r2_max_allowed = 0.4 * domain_diag;
if r2_avg > r2_max_allowed
    warning('r2_avg=%.3g exceeds domain size %.3g; capping to %.3g.', r2_avg, domain_diag, r2_max_allowed);
    r2_avg = r2_max_allowed;
end
% ---------- Window widths ----------
A   = (xhi-xlo)*(yhi-ylo);
rho_all = natom / max(A, epsr);

if useManual
    if useMixed
        dr1 = (2.355*sigr1); % +/- 1 FWHM
        dr2 = (2.355*sigr2);
    else
        %otherwise calculate the widths geometrically based on N widths
        dr1 = lam1*b*(2.355*sig1); % +/- 1 FWHM
        dr2 = lam2*b*(2.355*sig2);

    end
else
    % density-based (adaptive)
    dr1 = targetC1 / max(2*pi*rho_all*max(r1_avg,epsr), epsr);
    % For double network, size long window using sparse density to match larger spacing
    Nsparse_est = max(1, round(f_sparse * natom));
    % For double network, size long window using sparse density to match larger spacing
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

gap = 0.1*r1_avg;
r2_lower = r2_avg - dr2;
r2_upper = r2_avg + dr2;

if autoN2
    N2old = N2;
    R2AVG = 0.5*(r2_upper - r2_lower) + r2_lower;
    N2 = R2AVG/(lam2*b);
    fprintf("   Auto N2 mode enabled. Adjusting N2 from %.0d to %.0d\n",N2old,N2);
    options.bimodal.N2 = N2;
end

% ---------- Unpack atoms & maps ----------
ids = Atoms(:,1);
x   = Atoms(:,2);
y   = Atoms(:,3);

id2row = containers.Map('KeyType','int64','ValueType','int32');
for r=1:natom, id2row(int64(ids(r))) = int32(r); end

% ---------- Sparse subset (for double network long bonds) ----------
if double_network
    alpha = options.double_network.alpha; % desired spacing multiplier, e.g. 2.5
    avg_nn_spacing = sqrt((xhi-xlo)*(yhi-ylo)/natom); % crude estimate
    target_spacing = alpha * avg_nn_spacing;
    if alpha ==1
        target_spacing = 10;
    end
    isSparse = pick_uniform_sparse_nodes(x, y, f_sparse, target_spacing);
    sparse_idx = find(isSparse);
else
    isSparse = true(natom,1);
    sparse_idx = 1:natom;
end

% ---------- Grid (cell list) ----------
rmax = max(r1_upper, r2_upper);
if rmax<=0, rmax = max(xhi-xlo, yhi-ylo); end
cellSize = rmax;
nx = max(1, ceil((xhi-xlo)/cellSize));
ny = max(1, ceil((yhi-ylo)/cellSize));

cx = floor((x - xlo)/cellSize) + 1; cx = max(1, min(nx, cx));
cy = floor((y - ylo)/cellSize) + 1; cy = max(1, min(ny, cy));

Cells = cell(nx, ny);
for i=1:natom
    Cells{cx(i), cy(i)}(end+1) = i; %#ok<AGROW>
end

% ---------- Helpers ----------
% We track deg1 (type-1) and deg2 (type-2) separately
deg1 = zeros(natom,1);     % degrees from type-1 bonds
deg2 = zeros(natom,1);     % degrees from type-2 bonds

% store ROW indices first; convert to IDs at end
Btmp = zeros(Max_bond, 5);  % [bid, r1, r2, L, type]
nbond = 0;
countType2 = 0;

% adaptive per-pick multipliers
g_mul1 = 1.0; g_mul2 = 1.0;

% ---------- PASS ORDER ----------
tL = 0; tS = 0; tA = 0; tB = 0;

if long_first
    % ============================ PASS L: TYPE-2 FIRST ============================
    ntries = 0; no_progress = 0;
    tic
    while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
        if countType2 >= target_N2, break; end
        ntries = ntries + 1;

        % pick r1 (if double network, from sparse subset)
        if double_network
            r1 = sparse_idx(randi(numel(sparse_idx)));
        else
            r1 = randi(natom);
        end

        % capacity for long pass: keep deg2 < Max_peratom_bond (sanity)
        if deg2(r1) >= (Max_peratom_bond)
            no_progress = no_progress + 1; continue;
        end

        if double_network
            r2lo = r2_lower;
            r2hi = r2_upper;
        else
            dr2_pick = dr2 * g_mul2;
            % keep a dynamic gap vs short upper
            r2lo = max(r2_lower - (dr2 - dr2_pick), r1_upper + gap);
            r2hi = r2_upper + (dr2_pick - dr2);
        end
       
        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny,isPeriodic);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        % if double network: candidate partners must also be sparse
        if double_network
            neigh = neigh(isSparse(neigh));
            if isempty(neigh), no_progress = no_progress + 1; continue; end
        end

        neigh(neigh==r1) = [];
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        % partner must also have deg2 < Max_peratom_bond
        neigh = neigh(deg2(neigh) < (Max_peratom_bond));
        if isempty(neigh), no_progress = no_progress + 1; continue; end
        
        dxv = x(neigh) - x(r1);
        dyv = y(neigh) - y(r1);
        d = minimum_image(isPeriodic,dxv,dyv,Lx,Ly);

        in2 = (d >= r2lo) & (d <= r2hi);
        cand2 = neigh(in2);
        cand2 = exclude_existing_any(cand2, r1, Atoms, id2row); % exclude any existing edge (any type)

        if ~useManual
            % adaptive width
            C = numel(cand2);
            if C < Cmin
                g_mul2 = min(2.0, g_mul2*1.15);
                no_progress = no_progress + 1; continue;
            elseif C > Cmax
                g_mul2 = max(0.5, g_mul2*0.85);
            else
                g_mul2 = 1.0;
            end
        end

        if isempty(cand2), no_progress = no_progress + 1; continue; end

        % bias to lower deg2 partner
        [~,ord] = sort(deg2(cand2),'ascend');
        cand2 = cand2(ord);
        r2 = cand2( min(numel(cand2), randi( min(5, numel(cand2)) )) );

        L = minimum_image(isPeriodic,x(r2)-x(r1),y(r2)-y(r1),Lx,Ly);

        % add type-2 bond
        nbond = nbond + 1;
        Btmp(nbond,:) = [nbond, r1, r2, L, 2];
        countType2 = countType2 + 1;

        % Update per-type degrees
        deg2(r1) = deg2(r1) + 1;
        deg2(r2) = deg2(r2) + 1;

        % --- ensure neighbor list capacity before writing (avoid overflow) ---
        needCol_r1 = 5 + (Atoms(r1,5) + 1);
        needCol_r2 = 5 + (Atoms(r2,5) + 1);
        curCols = size(Atoms,2);
        if needCol_r1 > curCols || needCol_r2 > curCols
            growBy = max(16, max(needCol_r1, needCol_r2) - curCols);
            Atoms(:, curCols+1 : curCols+growBy) = 0;
        end
        % Update Atoms neighbor IDs for bookkeeping (IDs)
        Atoms(r1,5) = Atoms(r1,5) + 1;  Atoms(r1, 5+Atoms(r1,5)) = ids(r2);
        Atoms(r2,5) = Atoms(r2,5) + 1;  Atoms(r2, 5+Atoms(r2,5)) = ids(r1);

        no_progress = 0;
    end
    tL = toc;

    % ============================ PASS S: TYPE-1 AFTER ============================
    % Now place type-1 bonds as if type-2 didn't consume capacity:
    % i.e., enforce deg1 < Max_peratom_bond (ignore deg2).
    ntries = 0; no_progress = 0;
    tic
    while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
        % Optional: try to lift most nodes to >= z1_min on deg1
        if mean(deg1 >= z1_min) > 0.95
            % continue filling until bonds or stalls, but priority met
        end
        ntries = ntries + 1;

        r1 = randi(natom);
        if deg1(r1) >= Max_peratom_bond
            no_progress = no_progress + 1; continue;
        end

        % adaptive window for type-1
        dr1_pick = dr1 * g_mul1;
        rlo = max(r1_lower - (dr1 - dr1_pick), dmin + epsr);
        rhi = r1_upper + (dr1_pick - dr1);

        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        neigh(neigh==r1) = [];
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        % partner must have deg1 < Max_peratom_bond
        neigh = neigh(deg1(neigh) < Max_peratom_bond);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        dxv = x(neigh) - x(r1);
        dyv = y(neigh) - y(r1);
        d = minimum_image(isPeriodic,dxv,dyv,Lx,Ly);;

        in1 = (d >= rlo) & (d <= rhi);
        cand = neigh(in1);
        cand = exclude_existing_any(cand, r1, Atoms, id2row); % exclude any existing edge (long or short)

        if ~useManual
            % adaptive width
            C = numel(cand);
            if C < Cmin
                g_mul1 = min(2.0, g_mul1*1.15);
                no_progress = no_progress + 1; continue;
            elseif C > Cmax
                g_mul1 = max(0.5, g_mul1*0.85);
            else
                g_mul1 = 1.0;
            end
        end
        if isempty(cand), no_progress = no_progress + 1; continue; end

        r2 = cand(randi(numel(cand)));
        L = minimum_image(isPeriodic,x(r2)-x(r1),y(r2)-y(r1),Lx,Ly);

        % add type-1 bond
        nbond = nbond + 1;
        Btmp(nbond,:) = [nbond, r1, r2, L, 1];

        deg1(r1) = deg1(r1) + 1;
        deg1(r2) = deg1(r2) + 1;

        % --- ensure neighbor list capacity before writing ---
        needCol_r1 = 5 + (Atoms(r1,5) + 1);
        needCol_r2 = 5 + (Atoms(r2,5) + 1);
        curCols = size(Atoms,2);
        if needCol_r1 > curCols || needCol_r2 > curCols
            growBy = max(16, max(needCol_r1, needCol_r2) - curCols);
            Atoms(:, curCols+1 : curCols+growBy) = 0;
        end
        Atoms(r1,5) = Atoms(r1,5) + 1;  Atoms(r1, 5+Atoms(r1,5)) = ids(r2);
        Atoms(r2,5) = Atoms(r2,5) + 1;  Atoms(r2, 5+Atoms(r2,5)) = ids(r1);

        no_progress = 0;
    end
    tS = toc;

else
    % ---------- Original order (type-1 scaffold, then type-2) ----------
    % PASS A: type-1 scaffold (on total degree Atoms(:,5))
    ntries = 0; no_progress = 0;
    tic
    while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
        if mean(Atoms(:,5) >= z1_min) > 0.95, break; end
        ntries = ntries + 1;

        r1 = randi(natom);
        if Atoms(r1,5) >= Max_peratom_bond
            no_progress = no_progress + 1; continue;
        end

        dr1_pick = dr1 * g_mul1;
        rlo = max(r1_lower - (dr1 - dr1_pick), dmin + epsr);
        rhi = r1_upper + (dr1_pick - dr1);

        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        neigh(neigh==r1) = [];
        if isempty(neigh), no_progress = no_progress + 1; continue; end
        freeMask = (Atoms(neigh,5) < Max_peratom_bond);
        neigh = neigh(freeMask);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        dxv = x(neigh) - x(r1);
        dyv = y(neigh) - y(r1);
        d = minimum_image(isPeriodic,dxv,dyv,Lx,Ly);

        in1 = (d >= rlo) & (d <= rhi);
        cand = neigh(in1);

        cand = exclude_existing_any(cand, r1, Atoms, id2row);

        if ~useManual
            C = numel(cand);
            if C < Cmin
                g_mul1 = min(2.0, g_mul1*1.15); no_progress = no_progress + 1; continue;
            elseif C > Cmax
                g_mul1 = max(0.5, g_mul1*0.85);
            else
                g_mul1 = 1.0;
            end
        end
        if isempty(cand), no_progress = no_progress + 1; continue; end

        r2 = cand(randi(numel(cand)));
        L = minimum_image(isPeriodic,x(r2)-x(r1),y(r2)-y(r1),Lx,Ly);

        nbond = nbond + 1;
        Btmp(nbond,:) = [nbond, r1, r2, L, 1];

        % --- ensure neighbor list capacity before writing ---
        needCol_r1 = 5 + (Atoms(r1,5) + 1);
        needCol_r2 = 5 + (Atoms(r2,5) + 1);
        curCols = size(Atoms,2);
        if needCol_r1 > curCols || needCol_r2 > curCols
            growBy = max(16, max(needCol_r1, needCol_r2) - curCols);
            Atoms(:, curCols+1 : curCols+growBy) = 0;
        end
        Atoms(r1,5) = Atoms(r1,5) + 1;  Atoms(r1, 5+Atoms(r1,5)) = ids(r2);
        Atoms(r2,5) = Atoms(r2,5) + 1;  Atoms(r2, 5+Atoms(r2,5)) = ids(r1);

        no_progress = 0;
    end
    tA = toc;

    % PASS B: type-2
    ntries = 0; no_progress = 0;
    tic
    while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
        if (~useProb) && (countType2 >= target_N2), break; end
        if useProb && (countType2 > target_N2*1.1), break; end
        ntries = ntries + 1;

        % if double network, long bonds only among sparse nodes
        if double_network
            r1 = sparse_idx(randi(numel(sparse_idx)));
        else
            r1 = randi(natom);
        end
        if Atoms(r1,5) >= Max_peratom_bond
            no_progress = no_progress + 1; continue;
        end
        
        if double_network
            r2lo = r2_lower;
            r2hi = r2_upper;
        else
            dr2_pick = dr2 * g_mul2;
            % keep a dynamic gap vs short upper
            r2lo = max(r2_lower - (dr2 - dr2_pick), r1_upper + gap);
            r2hi = r2_upper + (dr2_pick - dr2);
        end

        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        if double_network
            neigh = neigh(isSparse(neigh));
            if isempty(neigh), no_progress = no_progress + 1; continue; end
        end

        neigh(neigh==r1) = [];
        if isempty(neigh), no_progress = no_progress + 1; continue; end
        freeMask = (Atoms(neigh,5) < Max_peratom_bond);
        neigh = neigh(freeMask);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        dxv = x(neigh) - x(r1);
        dyv = y(neigh) - y(r1);
        d = minimum_image(isPeriodic,dxv,dyv,Lx,Ly);

        in2 = (d >= r2lo) & (d <= r2hi);
        cand2 = neigh(in2);
        cand2 = exclude_existing_any(cand2, r1, Atoms, id2row);

        if ~useManual
            C = numel(cand2);
            if C < Cmin
                g_mul2 = min(2.0, g_mul2*1.15); no_progress = no_progress + 1; continue;
            elseif C > Cmax
                g_mul2 = max(0.5, g_mul2*0.85);
            else
                g_mul2 = 1.0;
            end
        end
        if isempty(cand2), no_progress = no_progress + 1; continue; end

        [~,ord] = sort(Atoms(cand2,5),'ascend');
        cand2 = cand2(ord);
        r2 = cand2( min(numel(cand2), randi( min(5, numel(cand2)) )) );

        L = minimum_image(isPeriodic,x(r2)-x(r1),y(r2)-y(r1),Lx,Ly);

        nbond = nbond + 1;
        Btmp(nbond,:) = [nbond, r1, r2, L, 2];
        countType2 = countType2 + 1;

        % --- ensure neighbor list capacity before writing ---
        needCol_r1 = 5 + (Atoms(r1,5) + 1);
        needCol_r2 = 5 + (Atoms(r2,5) + 1);
        curCols = size(Atoms,2);
        if needCol_r1 > curCols || needCol_r2 > curCols
            growBy = max(16, max(needCol_r1, needCol_r2) - curCols);
            Atoms(:, curCols+1 : curCols+growBy) = 0;
        end
        Atoms(r1,5) = Atoms(r1,5) + 1;  Atoms(r1, 5+Atoms(r1,5)) = ids(r2);
        Atoms(r2,5) = Atoms(r2,5) + 1;  Atoms(r2, 5+Atoms(r2,5)) = ids(r1);

        no_progress = 0;
    end
    tB = toc;
end

% ------------------ Trim and optional pruning ------------------
Btmp = Btmp(1:nbond,:);

if ~isempty(Btmp) && (min_keep > 0)
    changed = true;
    while changed
        % Recompute total degree from current Btmp
        deg_tot = zeros(natom,1);
        for k=1:size(Btmp,1)
            deg_tot(Btmp(k,2)) = deg_tot(Btmp(k,2)) + 1;
            deg_tot(Btmp(k,3)) = deg_tot(Btmp(k,3)) + 1;
        end
        to_del = find(deg_tot <= min_keep);
        if isempty(to_del), changed = false; break; end
        kill = false(size(Btmp,1),1);
        mark = false(natom,1); mark(to_del) = true;
        for k=1:size(Btmp,1)
            if mark(Btmp(k,2)) || mark(Btmp(k,3)), kill(k)=true; end
        end
        if any(kill), Btmp = Btmp(~kill,:); changed = true;
        else, changed = false;
        end
    end
end

% ------------------ Convert row indices to IDs; renumber ------------------
nb = size(Btmp,1);
BondsOut = zeros(nb,5);
for k=1:nb
    r1 = Btmp(k,2); r2 = Btmp(k,3);
    BondsOut(k,:) = [k, ids(r1), ids(r2), Btmp(k,4), Btmp(k,5)];
end

AtomsOut = Atoms;
Atoms    = AtomsOut;
Bonds    = BondsOut;

% ================== PRUNING: remove low-degree nodes + incident bonds ==================
% Requires: Atoms, Bonds are in ID-space (your AtomsOut/BondsOut); min_keep already set.
% R2016a-safe (no implicit expansion)

if min_keep > 0
    changed = true;
    while changed
        changed = false;

        if isempty(Bonds)
            % no bonds ? all degrees 0 < min_keep ? drop all atoms
            Atoms = []; Bonds = [];
            break
        end

        % Map atom IDs -> row indices
        id2row = containers.Map('KeyType','int64','ValueType','int32');
        for r = 1:size(Atoms,1)
            id2row(int64(Atoms(r,1))) = int32(r);
        end

        % Build degree from Bonds
        nA = size(Atoms,1);
        deg = zeros(nA,1);
        ri = zeros(size(Bonds,1),1,'int32');
        rj = ri;
        for k = 1:size(Bonds,1)
            ri(k) = id2row(int64(Bonds(k,2)));
            rj(k) = id2row(int64(Bonds(k,3)));
        end
        % count both endpoints
        idx   = [ri; rj];
        inc   = ones(numel(idx),1);
        deg   = deg + accumarray(double(idx), inc, [nA 1], @sum, 0);

        % nodes to delete this round
        delMask = (deg < min_keep);
        if ~any(delMask), break; end
        delIDs  = Atoms(delMask,1);

        % 1) delete bonds incident to any delIDs
        keepBond = ~ismember(Bonds(:,2), delIDs) & ~ismember(Bonds(:,3), delIDs);
        if any(~keepBond)
            Bonds = Bonds(keepBond,:);
            changed = true;
        end

        % 2) delete those atoms
        keepAtom = ~delMask;
        if any(~keepAtom)
            Atoms = Atoms(keepAtom,:);
            changed = true;
        end
    end
end

% If everything was pruned, stop here
if isempty(Atoms) || isempty(Bonds)
    % Ensure outputs are well-formed empty
    Atoms = zeros(0,5);   % [ID x y z num_bond], adjust if your layout differs
    Bonds = zeros(0,5);   % [bid i j L type]
    return
end

% ================== RENNUMBER IDs consecutively (LAMMPS-friendly) ==================
oldIDs      = Atoms(:,1);
newIDs      = (1:size(Atoms,1))';
Atoms(:,1)  = newIDs;

% map endpoints in Bonds from oldIDs -> newIDs
[tfI, locI] = ismember(Bonds(:,2), oldIDs);
[tfJ, locJ] = ismember(Bonds(:,3), oldIDs);
Bonds(tfI,2) = newIDs(locI(tfI));
Bonds(tfJ,3) = newIDs(locJ(tfJ));

% Reassign bond IDs to be consecutive 1..Nb (optional but clean)
if ~isempty(Bonds)
    Bonds(:,1) = (1:size(Bonds,1))';
end

% ================== REBUILD neighbor lists in Atoms(:,5:end) ==================
% wipe degree + neighbors
Atoms(:,5:end) = 0;

for k = 1:size(Bonds,1)
    ii = Bonds(k,2);
    jj = Bonds(k,3);
    % ii and jj are newIDs, which equal Atoms row indices after renumbering
    ri = ii; rj = jj;

    need_i = 5 + (Atoms(ri,5)+1);
    need_j = 5 + (Atoms(rj,5)+1);
    need   = max(need_i, need_j);

    % grow columns if needed (R2016a-safe)
    curC = size(Atoms,2);
    if need > curC
        Atoms(:, curC+1:need) = 0;
    end

    Atoms(ri,5) = Atoms(ri,5) + 1;  Atoms(ri, 5+Atoms(ri,5)) = jj;
    Atoms(rj,5) = Atoms(rj,5) + 1;  Atoms(rj, 5+Atoms(rj,5)) = ii;
end


% ------------------ Stats ------------------
tTot = tL + tS + tA + tB;
type1 = (Bonds(:,5)==1); type2 = (Bonds(:,5)==2);
fprintf('   LongFirst=%d  DoubleNet=%d (f_sparse=%.3f, Nsparse=%d) | T_long %.3fs, T_short %.3fs (Alt: T1 %.3fs, T2 %.3fs) | Total %.3fs\n', ...
        long_first, double_network, f_sparse, sum(isSparse), tL, tS, tA, tB, tTot);
fprintf('   Placed %d bonds with %d type1, %d type2\n', nb, sum(type1), sum(type2));

end

% ===== helpers =====
function neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic)
    Cx = cx(r1); Cy = cy(r1);
    neigh = [];
    
    if isPeriodic
        % gather neighbors for all adjacent cells periodic boundary conditions (with wrap)
        for dxCell=-1:1
            ix = Cx + dxCell;
            if ix < 1, ix = nx;
            elseif ix > nx, ix = 1;
            end
            for dyCell=-1:1
                iy = Cy + dyCell;
                if iy < 1, iy = ny;
                elseif iy > ny, iy = 1;
                end
                if ~isempty(Cells{ix,iy})
                    neigh = [neigh, Cells{ix,iy}]; %#ok<AGROW>
                end
            end
        end
    else
        % gather neighbors for all adjacent cells fixed boundary conditions (no wrap)
        for dxCell=-1:1
            ix = Cx + dxCell; if ix<1 || ix>nx, continue; end
            for dyCell=-1:1
                iy = Cy + dyCell; if iy<1 || iy>ny, continue; end
                if ~isempty(Cells{ix,iy})
                 neigh = [neigh, Cells{ix,iy}]; %#ok<AGROW>
                end
            end
        end
    end
end

function cand = exclude_existing_any(cand, r1, Atoms, id2row)
    % Exclude candidates already bonded to r1 (regardless of type).
    nb1 = Atoms(r1,5);
    if nb1 <= 0, return; end
    nbrIDs = Atoms(r1, 6:5+nb1);
    nbrIDs = nbrIDs(nbrIDs~=0);
    if isempty(nbrIDs), return; end
    nbrRows = zeros(size(nbrIDs));
    for jj=1:numel(nbrIDs)
        nbrRows(jj) = id2row(int64(nbrIDs(jj)));
    end
    if ~isempty(cand) && ~isempty(nbrRows)
        cand = setdiff(cand, nbrRows);
    end
end

function d = minimum_image(isPeriodic,dx,dy,Lx,Ly)
    
    %if fixed return Euclidean distance
    if ~isPeriodic
        d = sqrt(dx.^2 + dy.^2);
        return;
    end

    %otherwise apply minimum image convention
    dx_p = dx - Lx.*round(dx./Lx);
    dy_p = dy - Ly.*round(dy./Ly);
    d = sqrt(dx_p.^2 + dy_p.^2);
end
