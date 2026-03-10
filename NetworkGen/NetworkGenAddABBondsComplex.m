function [AtomsOut, BondsOut] = NetworkGenAddABBondsComplex(Domain, Atoms, Bonds, options)
% Add AB-only bonds (bond types 2 and 3): must connect atomType1 <-> atomType2.
% Enforces per-node cap on AB bonds: max_bondtype23_per_atom (total of type 2 + 3).
% Keeps existing bonds (type 1 backbone).

capAB = options.complex.max_bondtype23_per_atom;
fill  = options.complex.AB_fill;
fill  = max(0, min(1, fill));
f2    = options.complex.frac_bond2;
f2    = max(0, min(1, f2));

maxNbrTotal = Domain.maxNbrTotal;
typecol     = Domain.atomType_col;

natom = size(Atoms,1);
if natom == 0
    AtomsOut = Atoms;
    BondsOut = Bonds;
    return;
end

% Cutoff
if isfield(options,'complex') && isfield(options.complex,'Rcut_AB') && ~isempty(options.complex.Rcut_AB)
    Rcut = options.complex.Rcut_AB;
else
    Rcut = 1.85 * Domain.min_node_sep;
end
Rcut2 = Rcut*Rcut;

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;
Lx = xhi - xlo; Ly = yhi - ylo;

isPeriodic = strcmpi(options.boundary_box,'Periodic');

% Type sets
atype = Atoms(:,typecol);
is1 = (atype == 1);
is2 = (atype == 2);
N1 = sum(is1);
N2 = sum(is2);
if N1 == 0 || N2 == 0
    warning('No AB bonds added: need both atomType 1 and 2 present.');
    AtomsOut = Atoms;
    BondsOut = Bonds;
    return;
end

% Maximum possible AB bonds by capacity on each side
maxAB = min(capAB*N1, capAB*N2);
targetAB = round(fill * maxAB);
if targetAB <= 0
    AtomsOut = Atoms;
    BondsOut = Bonds;
    return;
end

target2 = round(f2 * targetAB);
target3 = targetAB - target2;

% Build adjacency from existing bonds to prevent duplicates
adj = spalloc(natom,natom,2*(size(Bonds,1)+targetAB+10));
for k=1:size(Bonds,1)
    i = Bonds(k,2); j = Bonds(k,3);
    adj(i,j) = 1; adj(j,i) = 1;
end

% Track AB degree per node (type 2 + 3 only)
degAB = zeros(natom,1);

% Build linked-cell grid
x = Atoms(:,2); y = Atoms(:,3);
hx = Rcut; hy = Rcut;
nx = max(1, floor((xhi-xlo)/hx));
ny = max(1, floor((yhi-ylo)/hy));
cx = floor((x - xlo)/hx) + 1; cx = max(1, min(nx, cx));
cy = floor((y - ylo)/hy) + 1; cy = max(1, min(ny, cy));
Cells = cell(nx, ny);
for i=1:natom
    Cells{cx(i), cy(i)}(end+1) = i; %#ok<AGROW>
end

% Preallocate added bonds
addedRows = zeros(targetAB, 5); % [bondID id1 id2 L btype] (bondID filled later)
nAdd = 0;

global_limit = max(200*targetAB, Domain.bond_global_try_limit);
stall_limit  = max(10*natom, Domain.max_attempts_without_progress);

ntries = 0; no_progress = 0;

while (nAdd < targetAB) && (ntries < global_limit) && (no_progress < stall_limit)
    ntries = ntries + 1;

    % choose bond type based on remaining counts
    rem2 = target2 - sum(addedRows(1:nAdd,5)==2);
    rem3 = target3 - sum(addedRows(1:nAdd,5)==3);
    if (rem2 + rem3) <= 0
        break;
    end
    if rem2 <= 0
        btype = 3;
    elseif rem3 <= 0
        btype = 2;
    else
        if rand < (rem2/(rem2+rem3))
            btype = 2;
        else
            btype = 3;
        end
    end

    r1 = randi(natom);
    if degAB(r1) >= capAB
        no_progress = no_progress + 1;
        continue;
    end

    % Must find opposite type
    t1 = atype(r1);
    if t1 ~= 1 && t1 ~= 2
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
    for kk=1:numel(neigh)
        r2 = neigh(kk);
        if r2 == r1, continue; end

        % Opposite type only
        if ~( (atype(r1)==1 && atype(r2)==2) || (atype(r1)==2 && atype(r2)==1) )
            continue;
        end

        if degAB(r2) >= capAB, continue; end
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

    nAdd = nAdd + 1;
    addedRows(nAdd,:) = [0, r1, r2, L, btype]; % bondID later

    adj(r1,r2) = 1; adj(r2,r1) = 1;
    degAB(r1) = degAB(r1) + 1;
    degAB(r2) = degAB(r2) + 1;

    no_progress = 0;
end

addedRows = addedRows(1:nAdd,:);
fprintf('   Added %d AB bonds (type2+3). Target was %d.\n', nAdd, targetAB);

% Append to Bonds, renumber bond IDs
BondsOut = Bonds;
if ~isempty(addedRows)
    nb0 = size(BondsOut,1);
    addedRows(:,1) = (nb0+1 : nb0+nAdd).';
    BondsOut = [BondsOut; addedRows];
end
BondsOut(:,1) = (1:size(BondsOut,1)).'; % ensure continuous IDs

% Rebuild neighbors in Atoms (all bonds)
AtomsOut = rebuild_neighbors_complex(Atoms, BondsOut, maxNbrTotal, typecol);

end

% ===== helpers =====
function At = rebuild_neighbors_complex(Atoms, Bonds, maxNbrTotal, typecol)
At = Atoms;
natom = size(At,1);

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