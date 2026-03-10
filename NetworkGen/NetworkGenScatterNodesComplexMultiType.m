function [Atoms] = NetworkGenScatterNodesComplexMultiType(Domain, options)
% Scatter nodes in 2D with minimum spacing, but also assign atomType=1/2 well-mixed.
% Output Atoms layout (compatible with legacy coords):
%   [ID x y z deg nbr1..nbrMax atomType]
%
% atomType stored at column Domain.atomType_col.

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;

Max_atom                  = Domain.Max_atom;
node_scatter_max_tries    = Domain.node_scatter_max_tries;
max_tries_per_node_sample = Domain.max_tries_per_node_sample;
min_node_sep2             = Domain.min_node_sep2;
dmin                      = sqrt(min_node_sep2);

% Total neighbor slots (backbone + AB)
if isfield(Domain,'maxNbrTotal')
    maxNbr = Domain.maxNbrTotal;
else
    % fallback
    maxNbr = Domain.Max_peratom_bond;
end
typecol = 6 + maxNbr;

% --------- Grid setup ---------
Lx = xhi - xlo;
Ly = yhi - ylo;
h  = dmin;
nx = max(1, ceil(Lx / h));
ny = max(1, ceil(Ly / h));

gridHeads = zeros(ny, nx, 'int32');
nextIdx   = zeros(Max_atom, 1, 'int32');

% --------- Storage ---------
Atoms  = zeros(Max_atom, typecol);
N_atom = 0;

% --------- Build shuffled type list for mixing ---------
phi2 = options.complex.phi_type2;
phi2 = max(0, min(1, phi2));
N2 = round(phi2 * Max_atom);
N1 = Max_atom - N2;
types = [ones(N1,1); 2*ones(N2,1)];
types = types(randperm(numel(types)));

    function [ci,cj] = coord2cell(x, y)
        cx = floor((x - xlo) / h) + 1;
        cy = floor((y - ylo) / h) + 1;
        if cx < 1, cx = 1; elseif cx > nx, cx = nx; end
        if cy < 1, cy = 1; elseif cy > ny, cy = ny; end
        ci = cx; cj = cy;
    end

    function ok = passes_minsep(xi, yi)
        [ci, cj] = coord2cell(xi, yi);
        ok = true;
        for dj = -1:1
            yj = cj + dj;
            if (yj < 1) || (yj > ny), continue; end
            for di = -1:1
                xi_cell = ci + di;
                if (xi_cell < 1) || (xi_cell > nx), continue; end
                head = gridHeads(yj, xi_cell);
                k = head;
                while k ~= 0
                    dx = xi - Atoms(k,2);
                    dy = yi - Atoms(k,3);
                    if (dx*dx + dy*dy) < min_node_sep2
                        ok = false;
                        return;
                    end
                    k = nextIdx(k);
                end
            end
        end
    end

    function insert_into_grid(idx)
        [ci, cj] = coord2cell(Atoms(idx,2), Atoms(idx,3));
        head = gridHeads(cj, ci);
        nextIdx(idx) = head;
        gridHeads(cj, ci) = int32(idx);
    end

% --------- Main loop ---------
global_scatter_tries = 0;
tic
while (N_atom < Max_atom) && (global_scatter_tries < node_scatter_max_tries)
    global_scatter_tries = global_scatter_tries + 1;

    accepted = false;
    per_node_tries = 0;

    while (~accepted) && (per_node_tries < max_tries_per_node_sample)
        per_node_tries = per_node_tries + 1;

        xi = xlo + Lx * rand;
        yi = ylo + Ly * rand;
        zi = 0;

        if N_atom == 0
            accepted = true;
        else
            accepted = passes_minsep(xi, yi);
        end
    end

    if ~accepted
        continue;
    end

    N_atom = N_atom + 1;
    Atoms(N_atom,1) = N_atom;
    Atoms(N_atom,2) = xi;
    Atoms(N_atom,3) = yi;
    Atoms(N_atom,4) = zi;
    Atoms(N_atom,5) = 0;              % degree
    Atoms(N_atom,6:5+maxNbr) = 0;     % neighbors
    Atoms(N_atom,typecol) = types(N_atom);

    insert_into_grid(N_atom);
end
fprintf('   Placed %d atoms (complex) in %4.4f sec \n', N_atom, toc);

if N_atom < Max_atom
    warning('Requested %d atoms, placed %d atoms.', Max_atom, N_atom);
    if N_atom == 0
        error('No atoms placed—aborting.');
    end
end

Atoms = Atoms(1:max(N_atom,1), :);

end