function Atoms_conn = NetworkBuildNeighborLists(Atoms, Bonds, max_nbr)
% -------------------------------------------------------------------------
% NetworkBuildNeighborLists
%
% Build per-node neighbor lists and degree from a Bonds list.
%
% INPUT:
%   Atoms : (Natoms x >=6)
%       col1 = id
%       col2 = x
%       col3 = y
%       col4 = z
%       col5 = type
%       col6 = isFixed  (from your lattice generator)
%
%   Bonds : (Nbonds x >=3)
%       col2 = atom i (id, not row)
%       col3 = atom j (id, not row)
%
%   max_nbr : maximum number of neighbors to store per node.
%             For your old style with 4 neighbors:
%                 max_nbr = 4
%             For a triangular lattice (degree up to ~6):
%                 max_nbr = 6
%
% OUTPUT:
%   Atoms_conn : [ID | X | Y | Z | num_bond | nbr1 | ... | nbr_max | spare]
%       'spare' here is isFixed from the original Atoms(:,6).
% -------------------------------------------------------------------------

if nargin < 3
    max_nbr = 4;   % default: match your old 4-neighbor layout
end

Natoms = size(Atoms,1);
Nbonds = size(Bonds,1);

% Pre-allocate degrees + neighbor lists (IDs)
deg  = zeros(Natoms,1);
nbrs = zeros(Natoms,max_nbr);

% Build id -> row index map (assumes IDs are 1..Natoms in order)
id_to_row = zeros(Natoms,1);
k = 1;
while k <= Natoms
    id_to_row(Atoms(k,1)) = k;
    k = k + 1;
end

% Loop over bonds and populate neighbor lists
b = 1;
while b <= Nbonds
    i_id = Bonds(b,2);  % atom id i
    j_id = Bonds(b,3);  % atom id j
    
    if i_id <= 0 || j_id <= 0
        b = b + 1;
        continue;
    end
    
    i_row = id_to_row(i_id);
    j_row = id_to_row(j_id);
    
    % Add j as neighbor of i
    di = deg(i_row);
    if di < max_nbr
        di = di + 1;
        nbrs(i_row, di) = j_id;
        deg(i_row) = di;
    end
    
    % Add i as neighbor of j
    dj = deg(j_row);
    if dj < max_nbr
        dj = dj + 1;
        nbrs(j_row, dj) = i_id;
        deg(j_row) = dj;
    end
    
    b = b + 1;
end

% Build output Atoms array:
% [ID | X | Y | Z | num_bond | nbr1..nbr_max | spare]
Atoms_conn = zeros(Natoms, 5 + max_nbr + 1);

% Copy positions
Atoms_conn(:,1) = Atoms(:,1);  % id
Atoms_conn(:,2) = Atoms(:,2);  % x
Atoms_conn(:,3) = Atoms(:,3);  % y
Atoms_conn(:,4) = Atoms(:,4);  % z

% num_bond
Atoms_conn(:,5) = deg;

% neighbor ids
Atoms_conn(:,6:(5+max_nbr)) = nbrs;

% spare: store isFixed (or anything else you like)
if size(Atoms,2) >= 6
    Atoms_conn(:,6+max_nbr) = Atoms(:,6);  % isFixed
else
    Atoms_conn(:,6+max_nbr) = 0;
end

end
