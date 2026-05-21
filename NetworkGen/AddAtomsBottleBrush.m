function [Atoms] = AddAtomsBottleBrush(obj)
% -------------------------------------------------------------------------
% AddAtomsBottleBrush
% - Scatter rod center-of-mass positions with minimum spacing
% - For each rod center, create a number of atoms in random orientation
% - Ensure rods don't intersect (periodic-aware)
%
% Rod length:
%   - obj.architecture.bottlebrush.rod_dispersity = false
%       every rod has exactly Nr atoms.
%   - obj.architecture.bottlebrush.rod_dispersity = true
%       each rod's atom count is drawn from the rod_length assignmentmode
%       (obj.architecture.bottlebrush.rod_length). See local_sample_rod_lengths
%       below for the supported distributions/methods.
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
    bb        = obj.architecture.bottlebrush;
    Nr        = bb.Nr;
    sigma_r   = bb.sigma_r;
    collision_multiplier = bb.rod_collision_multiplier;

    rod_dispersity = isfield(bb, 'rod_dispersity') && logical(bb.rod_dispersity);

    % --------- Boundary condition -----------------------------------------
    isPeriodic = strcmpi(obj.domain.boundary, 'periodic');

    % --------- Atoms column layout ----------------------------------------
    Max_peratom_bond = obj.peratom.Max_peratom_bond;
    ncols = 6 + Max_peratom_bond;  % ID | molID | X | Y | Z | deg | nbrs...

    % Calculate and log effective number of rods (junctions)
    max_possible_rods = floor(Max_atom / max(Nr,1));
    obj.log.print('   BottleBrush: Max_atom=%d, baseline Nr=%d per rod => max ~%d rods/junctions\n', ...
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

    obj.log.print('   BottleBrush: Max_atom=%d, baseline Nr=%d, max attempts=%d\n', ...
                  Max_atom, Nr, node_scatter_max_tries);

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

    % --------- Determine per-rod atom counts ------------------------------
    if rod_dispersity
        if ~isfield(bb, 'rod_length') || isempty(bb.rod_length)
            error(['AddAtomsBottleBrush: rod_dispersity is enabled but ' ...
                   'bottlebrush.rod_length (assignmentmode) is not set.']);
        end
        rod_lengths = local_sample_rod_lengths(bb.rod_length, N_rods);
        obj.log.print(['   Rod-length dispersity ON (mode=%s): ' ...
                       'mean=%.2f  min=%d  max=%d atoms/rod\n'], ...
                      lower(bb.rod_length.mode), mean(rod_lengths), ...
                      min(rod_lengths), max(rod_lengths));
    else
        rod_lengths = Nr * ones(N_rods, 1);
        obj.log.print('   Rod-length dispersity OFF: %d atoms/rod\n', Nr);
    end

    total_atoms = sum(rod_lengths);

    % --------- Expand rod centers to atoms ---------
    Atoms  = zeros(total_atoms, ncols);
    N_atom = 0;

    for i_rod = 1:N_rods

        Nr_i = rod_lengths(i_rod);

        x_center = RodCenters(i_rod, 1);
        y_center = RodCenters(i_rod, 2);
        z_center = RodCenters(i_rod, 3);

        dx = RodOrient(i_rod, 1);
        dy = RodOrient(i_rod, 2);
        dz = RodOrient(i_rod, 3);

        for i_atom = 1:Nr_i
            % Center the rod on its center-of-mass regardless of length
            s = (i_atom - 1 - (Nr_i - 1) / 2) * sigma_r;

            N_atom = N_atom + 1;

            Atoms(N_atom, 1) = N_atom;                % ID
            Atoms(N_atom, 2) = i_rod;                 % molID = rod index
            Atoms(N_atom, 3) = x_center + s * dx;     % X
            Atoms(N_atom, 4) = y_center + s * dy;     % Y
            Atoms(N_atom, 5) = z_center + s * dz;     % Z
            % col 6 = degree (0); cols 7..6+Max_peratom_bond = neighbors (0)
        end
    end

    obj.log.print('   Expanded to %d total atoms across %d rods\n', ...
                  N_atom, N_rods);

    Atoms = Atoms(1:N_atom, :);

end


% =========================================================================
% LOCAL: rod-length samplers
%
% These mirror the distribution conventions used by the AssignPerBond*
% family but produce per-rod integer atom counts (there is no bond length
% to couple to, so length-aligned methods are not applicable). Every rod is
% clamped to a minimum of MIN_ROD_LEN atoms so AddBondsRods always has at
% least one intra-rod bond to build.
% =========================================================================
function rod_lengths = local_sample_rod_lengths(am, N)

    MIN_ROD_LEN = 2;   % a rod needs >= 2 atoms to carry an internal bond

    rod_lengths = zeros(N, 1);
    if N == 0
        return;
    end

    mode = lower(am.mode);

    switch mode

        case 'mono'
            v = am.mono.value;
            if isempty(v)
                v = 20;
            end
            rod_lengths(:) = round(v);

        case 'uniform'
            lo = round(am.uniform.min_value);
            hi = round(am.uniform.max_value);
            if hi < lo
                tmp = lo; lo = hi; hi = tmp;
            end
            lo = max(lo, MIN_ROD_LEN);
            hi = max(hi, lo);
            rod_lengths = randi([lo, hi], N, 1);

        case {'poly', 'polydisperse'}
            rod_lengths = local_sample_poly(am.poly, N);

        case 'bimodal'
            rod_lengths = local_sample_bimodal(am.bimodal, N);

        otherwise
            error('AddAtomsBottleBrush: unknown rod_length.mode "%s".', am.mode);
    end

    rod_lengths = round(rod_lengths(:));
    rod_lengths = max(rod_lengths, MIN_ROD_LEN);

end


% -------------------------------------------------------------------------
% LOCAL: polydisperse rod-length sampling
%   Only 'pmf' and 'mono' methods are supported. 'geom' and 'range' are
%   bond-length-coupled and have no meaning for rod lengths.
% -------------------------------------------------------------------------
function vals = local_sample_poly(pd, N)

    method = lower(pd.method);

    min_N = 1;
    if isfield(pd, 'min_value') && ~isempty(pd.min_value)
        min_N = max(1, round(pd.min_value));
    end

    switch method

        case 'mono'
            if isfield(pd, 'pmf_mean') && ~isempty(pd.pmf_mean)
                v = round(pd.pmf_mean);
            else
                v = min_N;
            end
            vals = max(v, min_N) * ones(N, 1);

        case 'pmf'
            % Truncated geometric PMF on nu in [pmf_min, pmf_max] with a
            % target mean of pmf_mean. Same construction as the 'pmf' branch
            % of AssignPerBondPoly; the resulting value list is shuffled
            % across rods (there is no bond length to align to).
            nu0   = round(pd.pmf_min);
            nuMax = round(max(pd.pmf_min, pd.pmf_max));
            K     = nuMax - nu0;

            if K <= 0
                vals = max(nu0, min_N) * ones(N, 1);
                return;
            end

            targetMeanN = pd.pmf_mean;
            targetMeanK = max(0, targetMeanN - nu0);

            f = @(p) mean_k_of_p(p, K) - targetMeanK;

            p_lo = 1e-8;
            p_hi = 1 - 1e-8;

            Ps = linspace(1e-6, 1-1e-6, 200);
            Fs = zeros(size(Ps));
            for ii = 1:numel(Ps)
                Fs(ii) = f(Ps(ii));
            end

            bracket_found = false;
            a = NaN; btmp = NaN;
            for ii = 1:(numel(Ps)-1)
                if Fs(ii) == 0
                    a = Ps(ii); btmp = Ps(ii);
                    bracket_found = true;
                    break;
                elseif Fs(ii) * Fs(ii+1) < 0
                    a = Ps(ii); btmp = Ps(ii+1);
                    bracket_found = true;
                    break;
                end
            end

            if bracket_found
                if a == btmp
                    p_opt = a;
                else
                    p_opt = fzero(f, [a, btmp]);
                end
            else
                a_guess = max(eps, targetMeanK);
                p_opt = 1 / (a_guess + 1);
                p_opt = min(max(p_opt, p_lo), p_hi);
            end

            p = min(max(p_opt, p_lo), p_hi);
            r = 1 - p;

            denom = 1 - r^(K+1);
            Pk = (p * (r .^ (0:K)).') / denom;

            exp_counts  = Pk * N;
            base_counts = floor(exp_counts);
            remainder   = exp_counts - base_counts;

            assigned = sum(base_counts);
            deficit  = N - assigned;

            if deficit > 0
                [~, order] = sort(remainder, 'descend');
                for t = 1:deficit
                    base_counts(order(t)) = base_counts(order(t)) + 1;
                end
            elseif deficit < 0
                [~, order] = sort(remainder, 'ascend');
                for t = 1:(-deficit)
                    jj = order(t);
                    if base_counts(jj) > 0
                        base_counts(jj) = base_counts(jj) - 1;
                    end
                end
            end

            val_list = zeros(N, 1);
            ptr = 1;
            for k = 0:K
                cnt = base_counts(k+1);
                if cnt <= 0
                    continue;
                end
                val_list(ptr:ptr+cnt-1) = nu0 + k;
                ptr = ptr + cnt;
            end
            if ptr <= N
                val_list(ptr:N) = nu0 + K;
            end

            % No bond length to align to: shuffle across rods.
            vals = zeros(N, 1);
            rp = randperm(N).';
            vals(rp) = val_list;
            vals = max(vals, min_N);

        otherwise
            error(['AddAtomsBottleBrush: poly rod-length method "%s" is not ' ...
                   'supported for rod dispersity. Use "pmf" or "mono" ' ...
                   '(geom/range require bond lengths).'], pd.method);
    end

end


% -------------------------------------------------------------------------
% LOCAL: bimodal rod-length sampling
%   Each rod is assigned to population 1 or 2, then drawn from that
%   population. Only 'gaussian' and 'single' methods are supported.
% -------------------------------------------------------------------------
function vals = local_sample_bimodal(bd, N)

    method = lower(bd.method);

    min_N = 1;
    if isfield(bd, 'min_value') && ~isempty(bd.min_value)
        min_N = max(1, round(bd.min_value));
    end

    % Size of population 2 (e.g. the "long" population)
    if strcmpi(bd.height_mode, 'prob')
        n2 = round(bd.height_prob * N);
    else
        n2 = round(bd.height_count);
    end
    n2 = max(0, min(N, n2));
    n1 = N - n2;

    N1 = bd.mean_1;
    N2 = bd.mean_2;
    if isempty(N1), N1 = 35; end
    if isempty(N2), N2 = 60; end

    switch method

        case 'single'
            v1 = round(N1) * ones(n1, 1);
            v2 = round(N2) * ones(n2, 1);

        case 'gaussian'
            s1 = bd.std_1;
            s2 = bd.std_2;
            if isempty(s1) || s1 <= 0, s1 = max(1, 0.15 * N1); end
            if isempty(s2) || s2 <= 0, s2 = max(1, 0.15 * N2); end
            v1 = round(N1 + s1 * randn(n1, 1));
            v2 = round(N2 + s2 * randn(n2, 1));

        otherwise
            error(['AddAtomsBottleBrush: bimodal rod-length method "%s" is ' ...
                   'not supported for rod dispersity. Use "gaussian" or ' ...
                   '"single".'], bd.method);
    end

    vals = [v1; v2];
    vals = max(vals, min_N);

    % Shuffle so that rod index is uncorrelated with population
    vals = vals(randperm(N));
    vals = vals(:);

end


% -------------------------------------------------------------------------
% LOCAL: mean of a truncated geometric on k = 0..K with Pk ~ p (1-p)^k
% -------------------------------------------------------------------------
function mk = mean_k_of_p(p, K)

    r   = 1 - p;
    num = (1-p) .* (1 - (K+1)*r.^K + K*r.^(K+1)) ./ p;
    den = 1 - r.^(K+1);
    mk  = num ./ den;

end
