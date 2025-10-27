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

% --------- Unpack ---------
N_row             = size(Atoms,1);
Max_bond          = Domain.Max_bond;
Max_peratom_bond  = Domain.Max_peratom_bond;
global_limit      = Domain.bond_global_try_limit;
stall_limit       = Domain.max_attempts_without_progress;
min_keep          = Domain.min_degree_keep;

R_cut = Domain.min_node_sep * 2;
Rcut2 = Rcut*Rcut;

ids = Atoms(:,1);         % IDs by row
x   = Atoms(:,2);  y = Atoms(:,3);

% Build ID->row map (containers.Map is R2016a-safe)
id2row = containers.Map('KeyType','int64','ValueType','int32');
for r = 1:N_row
    id2row(int64(ids(r))) = int32(r);
end

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

    x1 = x(r1); y1 = y(r1);
    % brute-force candidate scan (R2016a robust)
    cand = zeros(16,1); ncan = 0;
    for r2 = 1:N_row
        if r2 == r1, continue; end
        if deg(r2) >= Max_peratom_bond, continue; end
        if adj(r1,r2) ~= 0, continue; end
        dx = x(r2)-x1; dy = y(r2)-y1;
        if (dx*dx + dy*dy) < Rcut2
            ncan = ncan + 1;
            if ncan > numel(cand), cand = [cand; zeros(numel(cand),1)]; end %#ok<AGROW>
            cand(ncan) = r2;
        end
    end

    if ncan == 0
        no_progress = no_progress + 1;
        continue;
    end

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
        dx = x(r2)-x(r1); dy = y(r2)-y(r1);
        BondsRows(k,3) = sqrt(dx*dx + dy*dy);
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

% Quick row lookup for writing neighbors (ids → row)
% (we already have id2row for robust mapping)

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
