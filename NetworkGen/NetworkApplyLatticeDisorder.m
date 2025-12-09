function Atoms = NetworkApplyLatticeDisorder(Atoms, Domain, options)
% -------------------------------------------------------------------------
% NetworkApplyLatticeDisorder
%
% Randomly perturbs the positions of NON-fixed atoms to introduce
% geometric disorder while keeping them inside the domain.
%
% disorder_level in [0,1] controls the strength:
%   0 -> no movement (perfect lattice)
%   1 -> max movement = disorder_max_frac_a * a
% -------------------------------------------------------------------------

% Extract lattice spacing
if ~isfield(options,'lattice') || ~isfield(options.lattice,'a')
    error('NetworkApplyLatticeDisorder: missing options.lattice.a');
end
a = options.lattice.a;

% Disorder level
disorder_level = 0;
if isfield(options.lattice,'disorder_level')
    disorder_level = options.lattice.disorder_level;
end
if disorder_level <= 0
    return;
end

% Max fraction of 'a' for displacement radius
disorder_max_frac_a = 0.35;
if isfield(options.lattice,'disorder_max_frac_a')
    disorder_max_frac_a = options.lattice.disorder_max_frac_a;
end

% Maximum displacement radius
r_max = disorder_level * disorder_max_frac_a * a;
if r_max <= 0
    return;
end

% Domain bounds
xlo = Domain.xlo;  xhi = Domain.xhi;
ylo = Domain.ylo;  yhi = Domain.yhi;

edgeTol = a*0.25;
if isfield(options,'lattice') && isfield(options.lattice,'edgeTol')
    edgeTol = options.lattice.edgeTol;
end

n = size(Atoms,1);

for i = 1:n
    % Skip fixed boundary nodes
    if Atoms(i,6) == 1
        continue;
    end
    
    % Sample random displacement uniformly in a disk of radius r_max
    % (sqrt(rand) to get uniform area distribution)
    rr    = r_max * sqrt(rand);
    theta = 2*pi*rand;
    dx    = rr * cos(theta);
    dy    = rr * sin(theta);
    
    x_new = Atoms(i,2) + dx;
    y_new = Atoms(i,3) + dy;
    
    % Clamp inside domain, leave a small edgeTol margin
    if x_new < xlo + edgeTol
        x_new = xlo + edgeTol;
    elseif x_new > xhi - edgeTol
        x_new = xhi - edgeTol;
    end
    
    if y_new < ylo + edgeTol
        y_new = ylo + edgeTol;
    elseif y_new > yhi - edgeTol
        y_new = yhi - edgeTol;
    end
    
    Atoms(i,2) = x_new;
    Atoms(i,3) = y_new;
end
end
