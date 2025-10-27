function [AtomsOut, BondsOut] = NetworkGenConnectNodesPolydisperse(Domain, Atoms, options)
% INPUT:
%   Atoms layout: [ ID | X | Y | Z | num_bond | nbr1 | nbr2 | nbr3 | nbr4 | ... ]
%   IDs are arbitrary (NOT equal to row index).
% OUTPUT:
%   BondsOut: [bondID | id1 | id2 | L]  (ids, not rows)
%   AtomsOut: Atoms with num_bond and neighbor slots rebuilt to match BondsOut
%
% Assumptions:
%   - Max_peratom_bond neighbor slots exist in Atoms (columns 6..(5+Max_peratom_bond))
%   - R2016a compatible, no implicit expansion
%   - Uses linked-cell (uniform grid) with cell size ~ Rcut

% --------- Unpack ---------
N_row             = size(Atoms,1);
Max_bond          = Domain.Max_bond;
Max_peratom_bond  = Domain.Max_peratom_bond;
global_limit      = Domain.bond_global_try_limit;
stall_limit       = Domain.max_attempts_without_progress;
min_keep          = Domain.min_degree_keep;

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;

% Cutoff selection (prefer explicit Domain.Rcut if provided)
if isfield(Domain,'Rcut')
    Rcut = Domain.Rcut;
else
    Rcut = 2 * Domain.min_node_sep;   % fallback
end
Rcut2 = Rcut*Rcut;

ids = Atoms(:,1);         % IDs by row
x   = Atoms(:,2);  y = Atoms(:,3);

% Build ID->row map (robust to id ~= row)
id2row = containers.Map('KeyType','int64','ValueType','int32');
for r = 1:N_row
    id2row(int64(ids(r))) = int32(r);
end

% --------- Build linked-cell grid (bins) ---------
% Cell size ~ Rcut; only search 3x3 neighboring cells per seed
hx = Rcut; hy = Rcut;
nx = max(1, floor((xhi - xlo)/hx));
ny = max(1, floor((yhi - ylo)/hy));

% Compute cell indices for each row
cx = floor((x - xlo) / hx) + 1;   % 1..nx
cy = floor((y - ylo) / hy) + 1;   % 1..ny
% Clamp to domain
cx(cx < 1) = 1; cx(cx > nx) = nx;
cy(cy < 1) = 1; cy(cy > ny) = ny;

% Linear bin index
binIdx = int32(cx + (cy-1)*nx);  % 1..nx*ny

% Build bins as cell array: bins{bin} = [row indices]
bins = accumarray(double(binIdx), (1:N_row)', [nx*ny, 1], @(v){v});

% --------- Bond creation in ROW space (store rows internally) ---------
deg = zeros(N_row,1);             % degrees by row
adj = sparse(N_row,N_row);        % 0/1 symmetric

BondsRows = zeros(Max_bond,3);    % [row1,row2,L]; bondID assigned later
nbond = 0;

no_progress = 0;
ntries = 0;

while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
    ntries = ntries + 1;

    % pick unsaturated row
    r1 = randi(N_row);
    if deg(r1) >= Max_peratom_bond
        no_progress = no_progress + 1;
        continue;
    end

    % 3x3 neighborhood bins around r1's cell
    c1x = cx(r1); c1y = cy(r1);

    % Gather candidate rows from neighbor cells
    candRows = []; % will grow; typical size small
    for dyc = -1:1
        yy = c1y + dyc;
        if (yy < 1) || (yy > ny), continue; end
        for dxc = -1:1
            xx = c1x + dxc;
            if (xx < 1) || (xx > nx), continue; end
            bId = xx + (yy-1)*nx;
            list = bins{bId};
            if ~isempty(list)
                candRows = [candRows; list]; %#ok<AGROW>
            end
        end
    end

    if isempty(candRows)
        no_progress = no_progress + 1;
        continue;
    end

    % Filter candidates: not self, unsaturated, not already connected, within Rcut
    x1 = x(r1); y1 = y(r1);
    cand = zeros(16,1); ncan = 0;

    % Iterate over local list only (fast)
    for idx = 1:numel(candRows)
        r2 = candRows(idx);
        if r2 == r1, continue; end
        if deg(r2) >= Max_peratom_bond, continue; end
        if adj(r1,r2) ~= 0, continue; end
        dxv = x(r2) - x1; dyv = y(r2) - y1;
        if (dxv*dxv + dyv*dyv) < Rcut2
            ncan = ncan + 1;
            if ncan > numel(cand), cand = [cand; zeros(numel(cand),1)]; end %#ok<AGROW>
            cand(ncan) = r2;
        end
    end

    if ncan == 0
        no_progress = no_progress + 1;
        continue;
    end

    % choose random neighbor among candidates
    r2 = cand(randi(ncan));
    L  = sqrt((x(r2)-x(r1))^2 + (y(r2)-y(r1))^2);

    nbond = nbond + 1;
    BondsRows(nbond,1) = r1;
    BondsRows(nbond,2) = r2;
    BondsRows(nbond,3) = L;

    deg(r1) = deg(r1) + 1; deg(r2) = deg(r2) + 1;
    adj(r1,r2) = 1; adj(r2,r1) = 1;

    no_progress = 0;
end

if ntries >= global_limit
    warning('Bond creation: hit global try limit (%d).', global_limit);
end
if no_progress >= stall_limit
    warning('Bond creation: local stall after %d attempts.', stall_limit);
end

BondsRows = BondsRows(1:nbond,:);

% --------- Iterative pruning (row space), no ID compaction ----------
if ~isempty(BondsRows)
    changed = true;
    while changed
        % recompute degree from current bonds
        deg(:) = 0;
        for k=1:size(BondsRows,1)
            deg(BondsRows(k,1)) = deg(BondsRows(k,1)) + 1;
            deg(BondsRows(k,2)) = deg(BondsRows(k,2)) + 1;
        end
        to_del = find(deg <= min_keep);
        if isempty(to_del)
            changed = false;
            break;
        end
        kill = false(size(BondsRows,1),1);
        mark = false(N_row,1); mark(to_del) = true;
        for k=1:size(BondsRows,1)
            if mark(BondsRows(k,1)) || mark(BondsRows(k,2))
                kill(k) = true;
            end
        end
        if any(kill)
            BondsRows = BondsRows(~kill,:);
            changed = true;
        else
            changed = false;
        end
    end

    % refresh lengths from coords (robust)
    for k=1:size(BondsRows,1)
        r1 = BondsRows(k,1); r2 = BondsRows(k,2);
        dxv = x(r2)-x(r1); dyv = y(r2)-y(r1);
        BondsRows(k,3) = sqrt(dxv*dxv + dyv*dyv);
    end
end

% --------- Build BondsOut in ID space; renumber bond IDs -----------
nb = size(BondsRows,1);
BondsOut = zeros(nb,4);
for k=1:nb
    r1 = BondsRows(k,1); r2 = BondsRows(k,2);
    id1 = ids(r1); id2 = ids(r2);
    L   = BondsRows(k,3);
    BondsOut(k,:) = [k, id1, id2, L];
end

% --------- Rebuild Atoms neighbors using **IDs** -------------------
% Clear neighbor metadata
Atoms(:,5) = 0;                             % num_bond
if size(Atoms,2) < 5+Max_peratom_bond
    % extend columns if needed (safeguard)
    Atoms(:, size(Atoms,2)+1 : 5+Max_peratom_bond) = 0;
else
    Atoms(:, 6 : 5+Max_peratom_bond) = 0;   % zero neighbor IDs
end

% Populate neighbor IDs from final bonds
for k=1:nb
    id1 = BondsOut(k,2);  id2 = BondsOut(k,3);
    r1  = id2row(int64(id1));
    r2  = id2row(int64(id2));

    % r1 side
    nb1 = Atoms(r1,5);
    if nb1 < Max_peratom_bond
        Atoms(r1,5) = nb1 + 1;
        Atoms(r1,5 + Atoms(r1,5)) = id2;   % store neighbor **ID**
    end
    % r2 side
    nb2 = Atoms(r2,5);
    if nb2 < Max_peratom_bond
        Atoms(r2,5) = nb2 + 1;
        Atoms(r2,5 + Atoms(r2,5)) = id1;   % store neighbor **ID**
    end
end

AtomsOut = Atoms;
end
