function [Atoms, latticeData] = NetworkGenLatticeScatterNodes(Domain, options)
% -------------------------------------------------------------------------
% NetworkGenLatticeScatterNodes
%
% Generate a perfect 2D hexagonal (triangular) lattice, then optionally
% apply geometric disorder to node positions.
% ...
% -------------------------------------------------------------------------

xlo = Domain.xlo;  xhi = Domain.xhi;
ylo = Domain.ylo;  yhi = Domain.yhi;

Lx = xhi - xlo;
Ly = yhi - ylo;

% Lattice spacing
if isfield(options,'lattice') && isfield(options.lattice,'a')
    a = options.lattice.a;
else
    error('NetworkGenLatticeScatterNodes: missing options.lattice.a');
end

edgeTol = a*0.25;
if isfield(options,'lattice') && isfield(options.lattice,'edgeTol')
    edgeTol = options.lattice.edgeTol;
end

dy =  a*sqrt(3)/2;

Ny_est = floor(Ly/dy) + 2;
Nx_est = floor(Lx/a) + 3;

idx_map = zeros(Ny_est, Nx_est);

maxNodes = Ny_est * Nx_est;
x_all = zeros(maxNodes,1);
y_all = zeros(maxNodes,1);

nat = 0;

for iy = 1:Ny_est
    y = ylo + (iy-1)*dy;
    if (y < ylo) || (y > yhi)
        continue;
    end
    
    % Even/odd row shift
    x_offset = 0;
    if mod(iy,2) == 1
        x_offset = a/2;
    end
    
    for ix = 1:Nx_est
        x = xlo + (ix-1)*a + x_offset;
        
        if (x < xlo) || (x > xhi)
            continue;
        end
        
        nat = nat + 1;
        x_all(nat) = x;
        y_all(nat) = y;
        idx_map(iy, ix) = nat;
    end
end

% Trim
x_all = x_all(1:nat);
y_all = y_all(1:nat);

% Fixed boundary nodes
isFixed = (x_all <= xlo + edgeTol) | (x_all >= xhi - edgeTol) | ...
          (y_all <= ylo + edgeTol) | (y_all >= yhi - edgeTol);

% Build Atoms matrix (perfect lattice)
Atoms = zeros(nat,6);
Atoms(:,1) = (1:nat).';   % id
Atoms(:,2) = x_all;       % x
Atoms(:,3) = y_all;       % y
Atoms(:,4) = 0;           % z
Atoms(:,5) = 1;           % type
Atoms(:,6) = double(isFixed);  % isFixed (0/1)

% ---- APPLY GEOMETRIC DISORDER (INTERIOR NODES ONLY) ----
if isfield(options,'lattice') && isfield(options.lattice,'disorder_level')
    if options.lattice.disorder_level > 0
        Atoms = NetworkApplyLatticeDisorder(Atoms, Domain, options);
    end
end

% Lattice info for connectivity
latticeData.idx_map = idx_map;
latticeData.Nx      = size(idx_map,2);
latticeData.Ny      = size(idx_map,1);
end
