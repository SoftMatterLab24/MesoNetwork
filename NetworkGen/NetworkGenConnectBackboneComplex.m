function [AtomsOut, BondsOut] = NetworkGenConnectBackboneComplex(Domain, Atoms, options)
% Build bondType=1 backbone, allowing any endpoint types (1-1,1-2,2-2),
% with per-atom backbone degree cap (default 5), then prune low-degree nodes.

natom = size(Atoms,1);

cap1  = options.complex.max_bondtype1_per_atom;
maxNbrTotal = Domain.maxNbrTotal;
typecol = Domain.atomType_col;

Max_bond_backbone = round(0.5*natom*cap1);

global_limit = Domain.bond_global_try_limit;
stall_limit  = Domain.max_attempts_without_progress;
min_keep     = Domain.min_degree_keep;

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;
Lx = xhi - xlo; Ly = yhi - ylo;

isPeriodic = strcmpi(options.boundary_box,'Periodic');

% Cutoff
if isfield(options,'complex') && isfield(options.complex,'Rcut_backbone') && ~isempty(options.complex.Rcut_backbone)
    Rcut = options.complex.Rcut_backbone;
else
    Rcut = 1.85 * Domain.min_node_sep;
end
Rcut2 = Rcut*Rcut;

x = Atoms(:,2); y = Atoms(:,3);

% --------- Build linked-cell grid ---------
hx = Rcut; hy = Rcut;
nx = max(1, floor((xhi-xlo)/hx));
ny = max(1, floor((yhi-ylo)/hy));

cx = floor((x - xlo)/hx) + 1; cx = max(1, min(nx, cx));
cy = floor((y - ylo)/hy) + 1; cy = max(1, min(ny, cy));

Cells = cell(nx, ny);
for i=1:natom
    Cells{cx(i), cy(i)}(end+1) = i; %#ok<AGROW>
end

% --------- Create bonds in row space ---------
deg1 = zeros(natom,1);
adj  = spalloc(natom,natom,2*max(1,Max_bond_backbone));

BondsRows = zeros(Max_bond_backbone,3); % [r1 r2 L]
nbond = 0; ntries = 0; no_progress = 0;

tic
while (nbond < Max_bond_backbone) && (no_progress < stall_limit) && (ntries < global_limit)
    ntries = ntries + 1;

    r1 = randi(natom);
    if deg1(r1) >= cap1
        no_progress = no_progress + 1;
        continue;
    end

    neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic);
    if isempty(neigh)
        no_progress = no_progress + 1;
        continue;
    end

    x1 = x(r1); y1 = y(r1);

    cand = zeros(32,1); ncan = 0;
    for kk = 1:numel(neigh)
        r2 = neigh(kk);
        if r2 == r1, continue; end
        if deg1(r2) >= cap1, continue; end
        if adj(r1,r2) ~= 0, continue; end

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
    L  = minimum_image(isPeriodic, x(r2)-x(r1), y(r2)-y(r1), Lx, Ly);

    nbond = nbond + 1;
    BondsRows(nbond,:) = [r1 r2 L];

    deg1(r1) = deg1(r1) + 1;
    deg1(r2) = deg1(r2) + 1;
    adj(r1,r2) = 1; adj(r2,r1) = 1;

    no_progress = 0;
end
fprintf('   Backbone: placed %d bonds in %4.4f sec \n', nbond, toc);

BondsRows = BondsRows(1:nbond,:);

% --------- Iterative pruning on backbone degrees ----------
pruned = false(natom,1);
if ~isempty(BondsRows)
    changed = true;
    while changed
        deg = zeros(natom,1);
        for k=1:size(BondsRows,1)
            deg(BondsRows(k,1)) = deg(BondsRows(k,1)) + 1;
            deg(BondsRows(k,2)) = deg(BondsRows(k,2)) + 1;
        end

        to_del = find(deg <= min_keep);
        if isempty(to_del)
            changed = false; break;
        end

        mark = false(natom,1); mark(to_del) = true;
        kill = false(size(BondsRows,1),1);
        for k=1:size(BondsRows,1)
            if mark(BondsRows(k,1)) || mark(BondsRows(k,2))
                kill(k) = true;
            end
        end

        if any(kill)
            BondsRows = BondsRows(~kill,:);
            pruned(to_del) = true;
            changed = true;
        else
            changed = false;
        end
    end
end

% Remove pruned atoms, remap rows
if any(pruned)
    keepMask = ~pruned;
    Atoms = Atoms(keepMask,:);
    old2new = zeros(numel(keepMask),1,'int32');
    old2new(keepMask) = int32(1:sum(keepMask));
    if ~isempty(BondsRows)
        keepBond = keepMask(BondsRows(:,1)) & keepMask(BondsRows(:,2));
        BondsRows = BondsRows(keepBond,:);
        BondsRows(:,1) = old2new(BondsRows(:,1));
        BondsRows(:,2) = old2new(BondsRows(:,2));
    end
end

natom = size(Atoms,1);
if natom == 0 || isempty(BondsRows)
    AtomsOut = zeros(0, 6+maxNbrTotal);
    BondsOut = zeros(0, 5);
    fprintf('   Backbone pruning removed everything.\n');
    return;
end

% Renumber IDs to match rows
Atoms(:,1) = (1:natom)';

% Build BondsOut with bondType column
nb = size(BondsRows,1);
BondsOut = zeros(nb,5);
for k=1:nb
    r1 = BondsRows(k,1);
    r2 = BondsRows(k,2);
    BondsOut(k,:) = [k, r1, r2, BondsRows(k,3), 1];
end

% Rebuild neighbor list (total slots, but currently only type1 bonds)
AtomsOut = rebuild_neighbors_complex(Atoms, BondsOut, maxNbrTotal, typecol);

fprintf('   Backbone final: %d atoms, %d bonds\n', size(AtomsOut,1), size(BondsOut,1));

end

% ===== helpers =====
function At = rebuild_neighbors_complex(Atoms, Bonds, maxNbrTotal, typecol)
At = Atoms;
natom = size(At,1);

% ensure columns exist
needCols = typecol;
if size(At,2) < needCols
    At(:, size(At,2)+1:needCols) = 0;
end

At(:,5) = 0;
At(:,6:5+maxNbrTotal) = 0;

for k=1:size(Bonds,1)
    i = Bonds(k,2);
    j = Bonds(k,3);
    di = At(i,5) + 1;
    if di <= maxNbrTotal
        At(i,5) = di;
        At(i,5+di) = j;
    end
    dj = At(j,5) + 1;
    if dj <= maxNbrTotal
        At(j,5) = dj;
        At(j,5+dj) = i;
    end
end
end

function neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny, isPeriodic)
Cx = cx(r1); Cy = cy(r1);
neigh = [];

if isPeriodic
    for dxCell=-1:1
        ix = Cx + dxCell;
        if ix < 1, ix = nx; elseif ix > nx, ix = 1; end
        for dyCell=-1:1
            iy = Cy + dyCell;
            if iy < 1, iy = ny; elseif iy > ny, iy = 1; end
            if ~isempty(Cells{ix,iy})
                neigh = [neigh, Cells{ix,iy}]; %#ok<AGROW>
            end
        end
    end
else
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

function d = minimum_image(isPeriodic, dx, dy, Lx, Ly)
if ~isPeriodic
    d = sqrt(dx*dx + dy*dy);
    return;
end
dxp = dx - Lx*round(dx/Lx);
dyp = dy - Ly*round(dy/Ly);
d = sqrt(dxp*dxp + dyp*dyp);
end