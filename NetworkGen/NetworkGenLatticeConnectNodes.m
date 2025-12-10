function [Atoms, Bonds] = NetworkGenLatticeConnectNodes(Domain, Atoms, latticeData, options)
% -------------------------------------------------------------------------
% NetworkGenLatticeConnectNodes
%
% Build first-neighbor bonds for a hexagonal (triangular) lattice.
%  - Uses the lattice idx_map (row,col -> atom id)
%  - No PBC
%  - L0 is ALWAYS the geometric distance between connected atoms,
%    so it is compatible with both perfect and disordered lattices.
%
% INPUT:
%   Domain (not really used, kept for consistency)
%   Atoms      : numeric matrix (Natoms x 6)
%                col1 = id
%                col2 = x
%                col3 = y
%                col4 = z
%                col5 = type
%                col6 = isFixed
%   latticeData.idx_map : Ny x Nx mapping (row,col) -> atom id (0 if none)
%
% OUTPUT:
%   Atoms : unchanged (passed through)
%   Bonds : numeric matrix (Nbonds x 5)
%           col1 = bond id
%           col2 = atom i
%           col3 = atom j
%           col4 = L0 (initial distance)
%           col5 = type (1)
% -------------------------------------------------------------------------

% Lattice index map
if ~isfield(latticeData,'idx_map')
    error('NetworkGenLatticeConnectNodes: latticeData.idx_map missing. Pass output from NetworkGenLatticeScatterNodes.');
end

idx_map = latticeData.idx_map;
Ny      = size(idx_map,1);
Nx      = size(idx_map,2);

Natoms  = size(Atoms,1);

% Estimated maximum number of bonds: ~3 bonds per node in this "forward" construction
maxBonds = 3 * Natoms;

i_list   = zeros(maxBonds,1);
j_list   = zeros(maxBonds,1);
L0_list  = zeros(maxBonds,1);
nb       = 0;

for iy = 1:Ny
    for ix = 1:Nx
        
        i_atom = idx_map(iy, ix);
        if i_atom == 0
            continue;
        end
        
        % 1) Horizontal neighbor: (iy, ix+1)
        if ix < Nx
            j_atom = idx_map(iy, ix+1);
            if j_atom > 0
                nb = nb + 1;
                i_list(nb) = i_atom;
                j_list(nb) = j_atom;
                
                dx = Atoms(j_atom,2) - Atoms(i_atom,2);  % x_j - x_i
                dy = Atoms(j_atom,3) - Atoms(i_atom,3);  % y_j - y_i
                L0_list(nb) = sqrt(dx*dx + dy*dy);
            end
        end
        
        % 2) Diagonal neighbors in row above (iy+1)
        if iy < Ny
            if mod(iy,2) == 1
                % Odd row: neighbors (iy+1, ix) and (iy+1, ix+1)
                j1 = idx_map(iy+1, ix);
                j2 = 0;
                if ix < Nx
                    j2 = idx_map(iy+1, ix+1);
                end
            else
                % Even row: neighbors (iy+1, ix) and (iy+1, ix-1)
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
                
                dx = Atoms(j1,2) - Atoms(i_atom,2);
                dy = Atoms(j1,3) - Atoms(i_atom,3);
                L0_list(nb) = sqrt(dx*dx + dy*dy);
            end
            
            if j2 > 0
                nb = nb + 1;
                i_list(nb) = i_atom;
                j_list(nb) = j2;
                
                dx = Atoms(j2,2) - Atoms(i_atom,2);
                dy = Atoms(j2,3) - Atoms(i_atom,3);
                L0_list(nb) = sqrt(dx*dx + dy*dy);
            end
        end
        
    end
end

% Trim arrays
i_list  = i_list(1:nb);
j_list  = j_list(1:nb);
L0_list = L0_list(1:nb);

% ---- Build numeric Bonds matrix compatible with your other code ----
% Use 5 columns:
%   1: bond id
%   2: atom i
%   3: atom j
%   4: L0
%   5: type
Bonds = zeros(nb,5);
Bonds(:,1) = (1:nb).';  % id
Bonds(:,2) = i_list;    % i
Bonds(:,3) = j_list;    % j
Bonds(:,4) = L0_list;   % reference length = geometric distance
Bonds(:,5) = 1;         % type (single bond type)

end
