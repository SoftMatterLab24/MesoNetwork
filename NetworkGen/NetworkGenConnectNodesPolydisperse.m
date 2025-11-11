function [AtomsOut, BondsOut] = NetworkGenConnectNodesPolydisperse(Domain, Atoms, options)
% NetworkGenConnectNodesPolydisperse - Connect Atoms nodes with Bonds
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
natom             = size(Atoms,1);
Max_bond          = Domain.Max_bond;
Max_peratom_bond  = Domain.Max_peratom_bond;
global_limit      = Domain.bond_global_try_limit;
stall_limit       = Domain.max_attempts_without_progress;
min_keep          = Domain.min_degree_keep;

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;
Lx = xhi - xlo; Ly = yhi - ylo;

% Mode options
isPeriodic  =  strcmpi(options.boundary_box,'Periodic');

% Cutoff selection (prefer explicit Domain.Rcut if provided)
if isfield(Domain,'Rcut')
    Rcut = Domain.Rcut;
else
    Rcut = 1.8 * Domain.min_node_sep;   % fallback
end
Rcut2 = Rcut*Rcut;

ids = Atoms(:,1);         % IDs by row
x   = Atoms(:,2);  y = Atoms(:,3);

% Build ID->row map (robust to id ~= row)
id2row = containers.Map('KeyType','int64','ValueType','int32');
for r = 1:natom
    id2row(int64(ids(r))) = int32(r);
end

% --------- Build linked-cell grid (bins) ---------
% Cell size ~ Rcut; only search 3x3 neighboring cells per seed
hx = Rcut; hy = Rcut;
nx = max(1, floor((xhi - xlo)/hx));
ny = max(1, floor((yhi - ylo)/hy));

cx = floor((x - xlo)/hx) + 1; cx = max(1, min(nx, cx));
cy = floor((y - ylo)/hy) + 1; cy = max(1, min(ny, cy));

Cells = cell(nx, ny);
for i=1:natom
    Cells{cx(i), cy(i)}(end+1) = i; %#ok<AGROW>
end
%% Compute cell indices for each row
%cx = floor((x - xlo) / hx) + 1;   % 1..nx
%cy = floor((y - ylo) / hy) + 1;   % 1..ny
%% Clamp to domain
%cx(cx < 1) = 1; cx(cx > nx) = nx;
%cy(cy < 1) = 1; cy(cy > ny) = ny;

%% Linear bin index
%binIdx = int32(cx + (cy-1)*nx);  % 1..nx*ny

%% Build bins as cell array: bins{bin} = [row indices]
%bins = accumarray(double(binIdx), (1:natom)', [nx*ny, 1], @(v){v});

% --------- Bond creation in ROW space (store rows internally) ---------
deg = zeros(natom,1);             % degrees by row
adj = sparse(natom,natom);        % 0/1 symmetric

BondsRows = zeros(Max_bond,3);    % [row1,row2,L]; bondID assigned later

nbond = 0;
ntries = 0; no_progress = 0;

tic
while (nbond < Max_bond) && (no_progress < stall_limit) && (ntries < global_limit)
    ntries = ntries + 1;

    % pick unsaturated row
    r1 = randi(natom);
    if deg(r1) >= Max_peratom_bond
        no_progress = no_progress + 1;
        continue;
    end

    %c1x = cx(r1); c1y = cy(r1);

    % Gather candidate rows from neighbor cells
    %candRows = []; % will grow; typical size small
    %for dyc = -1:1
    %    yy = c1y + dyc;
    %    if (yy < 1) || (yy > ny), continue; end
    %    for dxc = -1:1
    %        xx = c1x + dxc;
    %        if (xx < 1) || (xx > nx), continue; end
    %        bId = xx + (yy-1)*nx;
    %        list = bins{bId};
    %        if ~isempty(list)
    %            candRows = [candRows; list]; %#ok<AGROW>
    %        end
    %    end
    %end

    % 3x3 neighborhood bins around r1's cell
    neigh = gather_neighbors(r1, Cells, cx, cy, nx, ny,isPeriodic);
    if isempty(neigh), no_progress = no_progress + 1; continue; end

    % Filter candidates: not self, unsaturated, not already connected, within Rcut
    x1 = x(r1); y1 = y(r1);
    cand = zeros(16,1); ncan = 0;

    % Iterate over local list only (fast)
    for idx = 1:numel(neigh)
        r2 = neigh(idx);
        if r2 == r1, continue; end
        if deg(r2) >= Max_peratom_bond, continue; end
        if adj(r1,r2) ~= 0, continue; end
        dxv = x(r2) - x1; dyv = y(r2) - y1;
        d = minimum_image(isPeriodic,dxv,dyv,Lx,Ly);
        if (d*d) < Rcut2
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
    L = minimum_image(isPeriodic, x(r2)-x(r1), y(r2)-y(r1), Lx, Ly);
    
    nbond = nbond + 1;
    BondsRows(nbond,1) = r1;
    BondsRows(nbond,2) = r2;
    BondsRows(nbond,3) = L;

    deg(r1) = deg(r1) + 1; deg(r2) = deg(r2) + 1;
    adj(r1,r2) = 1; adj(r2,r1) = 1;

    no_progress = 0;
end
fprintf('   Placed %d bonds in %4.4f sec \n', nbond, toc);

if ntries >= global_limit
    warning('Bond creation: hit global try limit (%d).', global_limit);
end
if no_progress >= stall_limit
    warning('Bond creation: local stall after %d attempts.', stall_limit);
end

BondsRows = BondsRows(1:nbond,:);

% --------- Iterative pruning (remove low-degree nodes AND their bonds) ----------
pruned_atom_mask = false(natom,1);
if ~isempty(BondsRows)
    changed = true;
    while changed
        % recompute degree from current bonds (row space)
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

        % delete bonds incident to any low-degree node this round
        kill = false(size(BondsRows,1),1);
        mark = false(natom,1); mark(to_del) = true;
        for k = 1:size(BondsRows,1)
            if mark(BondsRows(k,1)) || mark(BondsRows(k,2)), kill(k) = true; end
        end
        if any(kill)
            BondsRows = BondsRows(~kill,:);
            pruned_atom_mask(to_del) = true;   % mark these atoms for removal
            changed = true;
        else
            changed = false;
        end
    end
end

% After pruning, drop pruned atoms from Atoms and compact rows
if any(pruned_atom_mask)
    Atoms = Atoms(~pruned_atom_mask, :);
    ids   = Atoms(:,1);          % update IDs vector to remaining atoms
    x     = Atoms(:,2);  y = Atoms(:,3);
    natom = size(Atoms,1);

    % Build old-row -> new-row remap for BondsRows endpoints
    old2new_row = zeros(numel(pruned_atom_mask),1,'int32');
    old2new_row(~pruned_atom_mask) = int32(1:natom);

    if ~isempty(BondsRows)
        % Remove any bonds that referenced deleted rows (safety)
        keepBond = ~pruned_atom_mask(BondsRows(:,1)) & ~pruned_atom_mask(BondsRows(:,2));
        BondsRows = BondsRows(keepBond,:);
        % Remap to new row indices
        BondsRows(:,1) = old2new_row(BondsRows(:,1));
        BondsRows(:,2) = old2new_row(BondsRows(:,2));
    end
end

% If everything got pruned, return empty well-formed outputs
if isempty(Atoms) || isempty(BondsRows)
    AtomsOut = zeros(0,5);       % [ID x y z num_bond]
    BondsOut = zeros(0,4);       % [bondID id1 id2 L]
    fprintf('   Pruned all atoms/bonds (min_keep=%d)\n', min_keep);
    return
end

% Refresh bond lengths from coordinates (robust, periodic-aware)
for k = 1:size(BondsRows,1)
    r1 = BondsRows(k,1); r2 = BondsRows(k,2);
    dxv = x(r2)-x(r1); dyv = y(r2)-y(r1);
    BondsRows(k,3) = minimum_image(isPeriodic, dxv, dyv, Lx, Ly);
end

% --------- Build BondsOut in ID space; renumber bond IDs -----------
nb = size(BondsRows,1);
BondsOut = zeros(nb,4);
for k = 1:nb
    r1 = BondsRows(k,1); r2 = BondsRows(k,2);
    id1 = ids(r1); id2 = ids(r2);
    L   = BondsRows(k,3);
    BondsOut(k,:) = [k, id1, id2, L];
end

% --------- Renumber atom IDs consecutively (clean for export) -------
oldIDs     = ids;
newIDs     = (1:natom)';          % consecutive
Atoms(:,1) = newIDs;

% Map bond endpoints oldID -> newID
[tfI,locI] = ismember(BondsOut(:,2), oldIDs);
[tfJ,locJ] = ismember(BondsOut(:,3), oldIDs);
BondsOut(tfI,2) = newIDs(locI(tfI));
BondsOut(tfJ,3) = newIDs(locJ(tfJ));

% Reassign bond IDs 1..Nb
if ~isempty(BondsOut), BondsOut(:,1) = (1:size(BondsOut,1))'; end

% --------- Rebuild Atoms neighbors using IDs ------------------------
% Ensure neighbor columns present
if size(Atoms,2) < 5+Max_peratom_bond
    Atoms(:, size(Atoms,2)+1 : 5+Max_peratom_bond) = 0;
end
% Clear degree + neighbor slots
Atoms(:,5) = 0;
Atoms(:,6:5+Max_peratom_bond) = 0;

% Since newIDs == row indices now, we can index directly
for k = 1:size(BondsOut,1)
    ii = BondsOut(k,2); jj = BondsOut(k,3);  % these equal row indices now
    % ii side
    nb1 = Atoms(ii,5) + 1;
    if nb1 <= Max_peratom_bond
        Atoms(ii,5) = nb1;
        Atoms(ii,5+nb1) = jj;   % store neighbor ID (which == row)
    end
    % jj side
    nb2 = Atoms(jj,5) + 1;
    if nb2 <= Max_peratom_bond
        Atoms(jj,5) = nb2;
        Atoms(jj,5+nb2) = ii;
    end
end

AtomsOut = Atoms;

fprintf('   Final: %d atoms, %d bonds (removed %d atoms, %d bonds by pruning; min_keep=%d)\n', ...
    size(AtomsOut,1), size(BondsOut,1), ...
    numel(pruned_atom_mask) - nnz(~pruned_atom_mask), ... % atoms removed
    nbond - size(BondsOut,1), min_keep);

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