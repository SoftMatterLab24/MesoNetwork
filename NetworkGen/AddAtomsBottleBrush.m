function [Atoms] = AddAtomsBottleBrush(obj)
% -------------------------------------------------------------------------
% AddAtomsBottleBrush
% - Scatter rod center-of-mass positions with minimum spacing
% - For each rod center, create Nr atoms in random orientation
% - Ensure rods don't intersect (periodic-aware)
%
% OUTPUT:
%   Atoms : [ID | molID | X | Y | Z | deg | nbr_1 ... nbr_{Max_peratom_bond}]
%           All atoms belonging to the same rod share the same molID (= rod index).
% -------------------------------------------------------------------------

    % --------- Unpack domain ---------
    xlo = obj.domain.xlo;
    xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;
    yhi = obj.domain.yhi;

    Max_atom                  = obj.domain.Max_atom;
    node_scatter_max_tries    = obj.domain.node_scatter_max_tries;
    max_tries_per_node_sample = obj.domain.max_tries_per_node_sample;
    min_node_sep2             = obj.domain.min_node_sep2;

    % --------- Unpack bottle-brush settings ---------
    Nr        = obj.architecture.bottlebrush.Nr;
    sigma_r   = obj.architecture.bottlebrush.sigma_r;
    collision_multiplier = obj.architecture.bottlebrush.rod_collision_multiplier;

    % --------- Boundary condition -----------------------------------------
    isPeriodic = strcmpi(obj.domain.boundary, 'periodic');

    % --------- Atoms column layout ----------------------------------------
    Max_peratom_bond = obj.peratom.Max_peratom_bond;
    ncols = 6 + Max_peratom_bond;  % ID | molID | X | Y | Z | deg | nbrs...

    % Calculate and log effective number of rods (junctions)
    max_possible_rods = floor(Max_atom / Nr);
    obj.log.print('   BottleBrush: Max_atom=%d, Nr=%d per rod => max ~%d rods/junctions\n', ...
                  Max_atom, Nr, max_possible_rods);

    dmin = sqrt(min_node_sep2);

    % Rod collision distance based on sigma_r (not rod length)
    rod_collision_dist = collision_multiplier * sigma_r;
    min_rod_sep2 = max(min_node_sep2, rod_collision_dist^2);

    % --------- Grid setup for rod center placement ---------
    Lx = xhi - xlo;
    Ly = yhi - ylo;

    h  = dmin;
    nx = max(1, ceil(Lx / h));
    ny = max(1, ceil(Ly / h));

    gridHeads = zeros(ny, nx, 'int32');
    nextIdx   = zeros(Max_atom, 1, 'int32');

    % Rod center positions (before expansion to atoms)
    RodCenters = zeros(Max_atom, 3);
    RodOrient  = zeros(Max_atom, 3);  % unit-length orientation vectors
    N_rods     = 0;

    % Final atom array (buffer room for one extra rod worth of atoms)
    Atoms  = zeros(Max_atom + Nr, ncols);
    N_atom = 0;

    % --------- Helper functions ---------
    function [ci, cj] = coord2cell(x, y)
        cx = floor((x - xlo) / h) + 1;
        cy = floor((y - ylo) / h) + 1;

        if cx < 1
            cx = 1;
        elseif cx > nx
            cx = nx;
        end

        if cy < 1
            cy = 1;
        elseif cy > ny
            cy = ny;
        end

        ci = cx;
        cj = cy;
    end

    function ok = passes_rod_minsep(xi, yi)
        [ci, cj] = coord2cell(xi, yi);
        ok = true;

        for dj = -1:1
            yj = cj + dj;
            if isPeriodic
                if yj < 1
                    yj = ny;
                elseif yj > ny
                    yj = 1;
                end
            else
                if (yj < 1) || (yj > ny)
                    continue;
                end
            end

            for di = -1:1
                xi_cell = ci + di;
                if isPeriodic
                    if xi_cell < 1
                        xi_cell = nx;
                    elseif xi_cell > nx
                        xi_cell = 1;
                    end
                else
                    if (xi_cell < 1) || (xi_cell > nx)
                        continue;
                    end
                end

                head = gridHeads(yj, xi_cell);
                k = head;

                while k ~= 0
                    dx = xi - RodCenters(k, 1);
                    dy = yi - RodCenters(k, 2);

                    if isPeriodic
                        dx = dx - Lx * round(dx / Lx);
                        dy = dy - Ly * round(dy / Ly);
                    end

                    if (dx*dx + dy*dy) < min_rod_sep2
                        ok = false;
                        return;
                    end

                    k = nextIdx(k);
                end
            end
        end
    end

    function insert_into_grid(idx)
        [ci, cj] = coord2cell(RodCenters(idx, 1), RodCenters(idx, 2));
        head = gridHeads(cj, ci);
        nextIdx(idx) = head;
        gridHeads(cj, ci) = int32(idx);
    end

    % --------- Place rod centers ---------
    global_scatter_tries = 0;
    tstart = tic;

    obj.log.print('   BottleBrush: Max_atom=%d, Nr=%d, max attempts=%d\n', Max_atom, Nr, node_scatter_max_tries);

    while (global_scatter_tries < node_scatter_max_tries)

        global_scatter_tries = global_scatter_tries + 1;

        accepted = false;
        per_node_tries = 0;

        while (~accepted) && (per_node_tries < max_tries_per_node_sample)

            per_node_tries = per_node_tries + 1;

            xi = xlo + Lx * rand;
            yi = ylo + Ly * rand;
            zi = 0.0;

            if N_rods == 0
                accepted = true;
            else
                accepted = passes_rod_minsep(xi, yi);
            end
        end

        if ~accepted
            continue;
        end

        N_rods = N_rods + 1;

        RodCenters(N_rods, 1) = xi;
        RodCenters(N_rods, 2) = yi;
        RodCenters(N_rods, 3) = zi;

        theta = 2 * pi * rand;
        RodOrient(N_rods, 1) = cos(theta);
        RodOrient(N_rods, 2) = sin(theta);
        RodOrient(N_rods, 3) = 0;

        insert_into_grid(N_rods);

    end

    obj.log.print('   Placed %d rod centers (%d junctions) in %4.4f sec (tries=%d)\n', ...
                  N_rods, N_rods, toc(tstart), global_scatter_tries);

    if N_rods < 1
        error('AddAtomsBottleBrush: No rod centers placed - aborting.');
    end

    % --------- Expand rod centers to atoms ---------
    for i_rod = 1:N_rods
        x_center = RodCenters(i_rod, 1);
        y_center = RodCenters(i_rod, 2);
        z_center = RodCenters(i_rod, 3);

        dx = RodOrient(i_rod, 1);
        dy = RodOrient(i_rod, 2);
        dz = RodOrient(i_rod, 3);

        for i_atom = 1:Nr
            s = (i_atom - 1 - (Nr - 1) / 2) * sigma_r;

            N_atom = N_atom + 1;

            Atoms(N_atom, 1) = N_atom;                % ID
            Atoms(N_atom, 2) = i_rod;                 % molID = rod index
            Atoms(N_atom, 3) = x_center + s * dx;     % X
            Atoms(N_atom, 4) = y_center + s * dy;     % Y
            Atoms(N_atom, 5) = z_center + s * dz;     % Z
            % col 6 = degree (0); cols 7..6+Max_peratom_bond = neighbors (0)
        end

        if N_atom >= 10*Max_atom
            obj.log.print('   Reached Max_atom=%d after placing rod %d (N_atom=%d)\n', Max_atom, i_rod, N_atom);
            break;
        end
    end

    obj.log.print('   Expanded to %d total atoms (%d rods x %d atoms/rod)\n', ...
                  N_atom, N_rods, Nr);

    Atoms = Atoms(1:N_atom, :);

end
