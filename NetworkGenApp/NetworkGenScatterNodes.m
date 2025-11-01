function [Atoms,obj] = NetworkGenScatterNodes(Domain,obj)
% NetworkGenScatterNodes - Scatter nodes in 2D with minimum spacing using a uniform grid (fast)
%
% INPUT:  Domain with fields:
%   xlo,xhi,ylo,yhi,zlo,zhi
%   Max_atom
%   node_scatter_max_tries
%   max_tries_per_node_sample
%   min_node_sep2   (squared minimum separation: d_min^2)
%
% OUTPUT:
%   Atoms: [ ID | X | Y | Z | num_bond | nbr1 | nbr2 | nbr3 | nbr4 | spare ]
%
% Notes:
% - R2016a compatible (no implicit expansion).
% - 2D scatter (Z=0).
% - Uses a spatial hash (uniform grid) to reduce neighbor checks.

% --------- Unpack ---------
xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;
% zlo = Domain.zlo; zhi = Domain.zhi;  % unused (planar)

Max_atom                  = Domain.Max_atom;
node_scatter_max_tries    = Domain.node_scatter_max_tries;
max_tries_per_node_sample = Domain.max_tries_per_node_sample;
min_node_sep2             = Domain.min_node_sep2;  % squared
dmin                      = sqrt(min_node_sep2);

% --------- Grid setup ---------
Lx = xhi - xlo;
Ly = yhi - ylo;
h  = dmin;                               % cell size ~ min separation
nx = max(1, ceil(Lx / h));
ny = max(1, ceil(Ly / h));

% Grid is stored as a "head index" per cell; 0 means empty.
% We keep a per-atom "next" array to form a linked list per cell.
gridHeads = zeros(ny, nx, 'int32');
nextIdx   = zeros(Max_atom, 1, 'int32');  % linked list next pointer

% --------- Storage for atoms ---------
Atoms  = zeros(Max_atom, 10);
N_atom = 0;

% --------- Helpers ---------
    function [ci,cj] = coord2cell(x, y)
        % map (x,y) to 1..nx,1..ny (clamped)
        cx = floor((x - xlo) / h) + 1;
        cy = floor((y - ylo) / h) + 1;
        if cx < 1, cx = 1; elseif cx > nx, cx = nx; end
        if cy < 1, cy = 1; elseif cy > ny, cy = ny; end
        ci = cx; cj = cy;
    end

    function ok = passes_minsep(xi, yi)
        % Check neighbors in 3x3 block around (ci,cj)
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
        % Insert atom idx into linked list of its cell
        [ci, cj]      = coord2cell(Atoms(idx,2), Atoms(idx,3));
        head          = gridHeads(cj, ci);
        nextIdx(idx)  = head;
        gridHeads(cj, ci) = int32(idx);
    end

% --------- Main loop ---------
global_scatter_tries = 0;
tic
while (N_atom < Max_atom) && (global_scatter_tries < node_scatter_max_tries)
    global_scatter_tries = global_scatter_tries + 1;

    % Try to place ONE point this outer iteration
    accepted = false;
    per_node_tries = 0;

    while (~accepted) && (per_node_tries < max_tries_per_node_sample)
        per_node_tries = per_node_tries + 1;

        xi = xlo + Lx * rand;
        yi = ylo + Ly * rand;
        zi = 0;

        if N_atom == 0
            accepted = true;  % first always accepted
        else
            accepted = passes_minsep(xi, yi);
        end
    end

    if ~accepted
        continue;  % try again with a new random candidate
    end

    % Place atom
    N_atom = N_atom + 1;
    if mod(N_atom, 1000) == 0
        N_atom; %#ok<NOPRT>
    end
    Atoms(N_atom,1) = N_atom;
    Atoms(N_atom,2) = xi;
    Atoms(N_atom,3) = yi;
    Atoms(N_atom,4) = zi;

    % Insert into grid
    insert_into_grid(N_atom);
end
newline = sprintf('   Placed %d atoms in %4.4f sec \n', N_atom, toc);
obj.log = append(obj.log, newline);

fprintf('   Placed %d atoms in %4.4f sec \n', N_atom, toc);

% Trim / warn
if N_atom < Max_atom
    newline = sprintf('   Warning: Requested %d atoms, placed %d atoms.\n', Max_atom, N_atom);
    obj.log = append(obj.log, newline);
    warning('Requested %d atoms, placed %d atoms.', Max_atom, N_atom);
    if N_atom == 0
        newline = sprintf('Error: No atoms placed—aborting.\n');
        obj.log = append(obj.log, newline);
        error('No atoms placed—aborting.');
    end
end
Atoms = Atoms(1:max(N_atom,1), :);

end
