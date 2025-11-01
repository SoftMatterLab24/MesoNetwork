function [Atoms,Bonds] = NetworkGenConnectNodesBimodal(Domain,Atoms,options,advancedOptions)
% NetworkGenConnectNodesBimodal (adaptive, grid-accelerated)
% Adds 'long_first' mode: place type-2 bonds first (fixed N2_number),
% then place type-1 bonds as if type-2 bonds don't consume capacity.
%
% Returns Bonds(:,5) = 1 or 2 for type.

% ---------- Unpack domain & defaults ----------
natom            = size(Atoms,1);
Max_bond         = Domain.Max_bond;
Max_peratom_bond = Domain.Max_peratom_bond;
xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;

if isfield(Domain,'bond_global_try_limit'), global_limit = Domain.bond_global_try_limit; else, global_limit = 5e7; end
if isfield(Domain,'max_attempts_without_progress'), stall_limit = Domain.max_attempts_without_progress; else, stall_limit = 5e5; end
if isfield(Domain,'min_degree_keep'), min_keep = Domain.min_degree_keep; else, min_keep = 0; end
if isfield(Domain,'min_node_sep'), dmin = Domain.min_node_sep; else, dmin = 0; end
epsr = 1e-9;

% ---------- Options ----------
b   = options.b;
N1  = options.bimodal.N1;
N2  = options.bimodal.N2;

useProb     = strcmpi(options.bimodal.distribution_height_mode,'prob');
long_first  = isfield(options.bimodal,'long_first') && options.bimodal.long_first;

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

% ---------- Refine r1_avg and r2_avg based on geometry ----------
r1_avg = b*sqrt(N1);
r2_avg = b*sqrt(N2);

r_min_allowed = max(dmin * 1.8, b * 0.5);
if r1_avg < r_min_allowed
    warning('r1_avg=%.3g below min separation %.3g; lifting to %.3g.', r1_avg, dmin, r_min_allowed);
    r1_avg = r_min_allowed;
end
if r2_avg < 1.8 * r1_avg
    warning('r2_avg too close to r1_avg (%.3g vs %.3g). Adjusting r2_avg.', r2_avg, r1_avg);
    r2_avg = 1.8 * r1_avg;
end

Lx = xhi - xlo; Ly = yhi - ylo;
domain_diag = hypot(Lx, Ly);
r2_max_allowed = 0.4 * domain_diag;
if r2_avg > r2_max_allowed
    warning('r2_avg=%.3g exceeds domain size %.3g; capping to %.3g.', r2_avg, domain_diag, r2_max_allowed);
    r2_avg = r2_max_allowed;
end

% ---------- Density-based window widths ----------
A   = (xhi-xlo)*(yhi-ylo);
rho = natom / max(A, eps);

dr1 = targetC1 / max(2*pi*rho*max(r1_avg,eps), eps);
dr2 = targetC2 / max(2*pi*rho*max(r2_avg,eps), eps);

r1_lower = max(r1_avg - dr1, dmin + epsr);
r1_upper = r1_avg + dr1;

gap = 0.1*r1_avg;
r2_lower = max(r2_avg - dr2, r1_upper + gap);
r2_upper = r2_avg + dr2;

% ---------- Unpack atoms & maps ----------
ids = Atoms(:,1);
x   = Atoms(:,2);
y   = Atoms(:,3);

id2row = containers.Map('KeyType','int64','ValueType','int32');
for r=1:natom, id2row(int64(ids(r))) = int32(r); end

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

        % capacity for long pass: keep deg2 < Max_peratom_bond (sanity)
        r1 = randi(natom);
        if deg2(r1) >= (Max_peratom_bond-1)
            no_progress = no_progress + 1; continue;
        end

        dr2_pick = dr2 * g_mul2;
        r2lo = max(r2_lower - (dr2 - dr2_pick), r1_upper + 0);
        r2hi = r2_upper + (dr2_pick - dr2);

        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        neigh(neigh==r1) = [];
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        % partner must also have deg2 < Max_peratom_bond
        neigh = neigh(deg2(neigh) < (Max_peratom_bond-1));
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        dxv = x(neigh) - x(r1);
        dyv = y(neigh) - y(r1);
        d   = sqrt(dxv.^2 + dyv.^2);

        in2 = (d >= r2lo) & (d <= r2hi);
        cand2 = neigh(in2);
        cand2 = exclude_existing_any(cand2, r1, Atoms, id2row); % exclude any existing edge (any type)

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
        if isempty(cand2), no_progress = no_progress + 1; continue; end

        % bias to lower deg2 partner
        [~,ord] = sort(deg2(cand2),'ascend');
        cand2 = cand2(ord);
        r2 = cand2( min(numel(cand2), randi( min(5, numel(cand2)) )) );

        L  = hypot(x(r2)-x(r1), y(r2)-y(r1));

        % add type-2 bond
        nbond = nbond + 1;
        Btmp(nbond,:) = [nbond, r1, r2, L, 2];
        countType2 = countType2 + 1;

        % Update per-type degrees
        deg2(r1) = deg2(r1) + 1;
        deg2(r2) = deg2(r2) + 1;

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

        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        neigh(neigh==r1) = [];
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        % partner must have deg1 < Max_peratom_bond
        neigh = neigh(deg1(neigh) < Max_peratom_bond);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        dxv = x(neigh) - x(r1);
        dyv = y(neigh) - y(r1);
        d   = sqrt(dxv.^2 + dyv.^2);

        in1 = (d >= rlo) & (d <= rhi);
        cand = neigh(in1);
        cand = exclude_existing_any(cand, r1, Atoms, id2row); % exclude any existing edge (long or short)

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
        if isempty(cand), no_progress = no_progress + 1; continue; end

        r2 = cand(randi(numel(cand)));
        L  = hypot(x(r2)-x(r1), y(r2)-y(r1));

        % add type-1 bond
        nbond = nbond + 1;
        Btmp(nbond,:) = [nbond, r1, r2, L, 1];

        deg1(r1) = deg1(r1) + 1;
        deg1(r2) = deg1(r2) + 1;

        Atoms(r1,5) = Atoms(r1,5) + 1;  Atoms(r1, 5+Atoms(r1,5)) = ids(r2);
        Atoms(r2,5) = Atoms(r2,5) + 1;  Atoms(r2, 5+Atoms(r2,5)) = ids(r1);

        no_progress = 0;
    end
    tS = toc;

else
    % ---------- Your original order (type-1 scaffold, then type-2) ----------
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

        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        neigh(neigh==r1) = [];
        if isempty(neigh), no_progress = no_progress + 1; continue; end
        freeMask = (Atoms(neigh,5) < Max_peratom_bond);
        neigh = neigh(freeMask);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        dxv = x(neigh) - x(r1);
        dyv = y(neigh) - y(r1);
        d   = sqrt(dxv.^2 + dyv.^2);

        in1 = (d >= rlo) & (d <= rhi);
        cand = neigh(in1);

        cand = exclude_existing_any(cand, r1, Atoms, id2row);

        C = numel(cand);
        if C < Cmin
            g_mul1 = min(2.0, g_mul1*1.15); no_progress = no_progress + 1; continue;
        elseif C > Cmax
            g_mul1 = max(0.5, g_mul1*0.85);
        else
            g_mul1 = 1.0;
        end
        if isempty(cand), no_progress = no_progress + 1; continue; end

        r2 = cand(randi(numel(cand)));
        L  = hypot(x(r2)-x(r1), y(r2)-y(r1));

        nbond = nbond + 1;
        Btmp(nbond,:) = [nbond, r1, r2, L, 1];

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

        r1 = randi(natom);
        if Atoms(r1,5) >= Max_peratom_bond
            no_progress = no_progress + 1; continue;
        end

        dr2_pick = dr2 * g_mul2;
        r2lo = max(r2_lower - (dr2 - dr2_pick), r1_upper + 0);
        r2hi = r2_upper + (dr2_pick - dr2);

        neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        neigh(neigh==r1) = [];
        if isempty(neigh), no_progress = no_progress + 1; continue; end
        freeMask = (Atoms(neigh,5) < Max_peratom_bond);
        neigh = neigh(freeMask);
        if isempty(neigh), no_progress = no_progress + 1; continue; end

        dxv = x(neigh) - x(r1);
        dyv = y(neigh) - y(r1);
        d   = sqrt(dxv.^2 + dyv.^2);

        in2 = (d >= r2lo) & (d <= r2hi);
        cand2 = neigh(in2);
        cand2 = exclude_existing_any(cand2, r1, Atoms, id2row);

        C = numel(cand2);
        if C < Cmin
            g_mul2 = min(2.0, g_mul2*1.15); no_progress = no_progress + 1; continue;
        elseif C > Cmax
            g_mul2 = max(0.5, g_mul2*0.85);
        else
            g_mul2 = 1.0;
        end
        if isempty(cand2), no_progress = no_progress + 1; continue; end

        [~,ord] = sort(Atoms(cand2,5),'ascend');
        cand2 = cand2(ord);
        r2 = cand2( min(numel(cand2), randi( min(5, numel(cand2)) )) );

        L  = hypot(x(r2)-x(r1), y(r2)-y(r1));

        nbond = nbond + 1;
        Btmp(nbond,:) = [nbond, r1, r2, L, 2];
        countType2 = countType2 + 1;

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

% ------------------ Rebuild Atoms neighbor lists from BondsOut -----------
% Atoms(:,5) = 0;
% if size(Atoms,2) < 5+Max_peratom_bond
%     Atoms(:, size(Atoms,2)+1 : 5+Max_peratom_bond) = 0;
% else
%     Atoms(:, 6 : 5+Max_peratom_bond) = 0;
% end
% for k=1:nb
%     id1 = BondsOut(k,2); id2 = BondsOut(k,3);
%     r1  = id2row(int64(id1)); r2 = id2row(int64(id2));
%     nb1 = Atoms(r1,5);
% %     if nb1 < Max_peratom_bond
%         Atoms(r1,5) = nb1 + 1; Atoms(r1,5+Atoms(r1,5)) = id2;
% %     end
%     nb2 = Atoms(r2,5);
% %     if nb2 < Max_peratom_bond
%         Atoms(r2,5) = nb2 + 1; Atoms(r2,5+Atoms(r2,5)) = id1;
% %     end
% end

AtomsOut = Atoms;
Atoms    = AtomsOut;
Bonds    = BondsOut;

% ------------------ Stats ------------------
tTot = tL + tS + tA + tB;
type1 = (Bonds(:,5)==1); type2 = (Bonds(:,5)==2);
fprintf('   LongFirst=%d | T_long %.3fs, T_short %.3fs (Alt: T1 %.3fs, T2 %.3fs) | Total %.3fs\n', ...
        long_first, tL, tS, tA, tB, tTot);
fprintf('   Bonds: %d (type1=%d, type2=%d)\n', nb, sum(type1), sum(type2));

end

% ===== helpers =====
function neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny)
    Cx = cx(r1); Cy = cy(r1);
    neigh = [];
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
