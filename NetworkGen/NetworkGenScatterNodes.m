function [Atoms] = NetworkGenScatterNodes(Domain)
% NetworkGenScatterNodes - Scatter nodes in the domain with minimum spacing
%
% INPUT:
%   Domain: structure with fields
%   
% OUTPUT:
%   Atoms: [ ID | X | Y | Z | num_bond | nbr1 | nbr2 | nbr3 | nbr4 | spare ]
%
% Assumptions:
%   - R2016a compatible, no implicit expansion
%   - Scatter in 2D plane (Z=0)


% --------- Unpack ---------
xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;
zlo = Domain.zlo; zhi = Domain.zhi;

Max_atom                        = Domain.Max_atom;
node_scatter_max_tries          = Domain.node_scatter_max_tries;
max_tries_per_node_sample       = Domain.max_tries_per_node_sample;
min_node_sep2                   = Domain.min_node_sep2;
Rmax                            = min_node_sep2;

% --------- Setup ---------
Atoms  = zeros(Max_atom, 10);
N_atom = 0;

global_scatter_tries = 0;
tic
while (N_atom < Max_atom) && (global_scatter_tries < node_scatter_max_tries)
    global_scatter_tries = global_scatter_tries + 1;

    % Try to draw a point that respects min_node_sep from all existing nodes
    accepted = false;
    per_node_tries = 0;

    while (~accepted) && (per_node_tries < max_tries_per_node_sample)
        per_node_tries = per_node_tries + 1;

        xi = xlo + (xhi - xlo) * rand;
        yi = ylo + (yhi - ylo) * rand;
        zi = 0; % planar

        if N_atom == 0
            accepted = true;  % first node always accepted
        else
            % Check squared distances to all currently placed nodes
            ok = true;
            for k = 1:N_atom
                dx = xi - Atoms(k,2);
                dy = yi - Atoms(k,3);
                if (dx*dx + dy*dy) < min_node_sep2
                    ok = false;
                    break;
                end
            end
            accepted = ok;
        end
    end

    if ~accepted
        % Could not find a valid spot for this node within per-node tries;
        % move on (global loop will stop when global_scatter_tries hits cap).
        continue;
    end

    % Place node
    N_atom = N_atom + 1;
    Atoms(N_atom,1) = N_atom;
    Atoms(N_atom,2) = xi;
    Atoms(N_atom,3) = yi;
    Atoms(N_atom,4) = zi;
end
fprintf('   Placed %d atoms in %4.4f sec \n', N_atom, toc);

% Warning and error checks
if N_atom < Max_atom
    warning('Requested %d atoms, placed %d atoms.', ...
            Max_atom, N_atom);
    if N_atom == 0
        error('No atoms placed—aborting.');
    end
    Atoms = Atoms(1:N_atom, :);
end

end