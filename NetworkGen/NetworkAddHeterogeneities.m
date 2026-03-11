function [AtomsOut, BondsOut] = NetworkAddHeterogeneities(Atoms, Bonds, Domain, options)
% NetworkAddHeterogeneities  Insert heterogeneous void defects into a network
%
% Removes crosslink nodes (atoms) and their bonds that fall inside randomly
% placed void regions, mimicking sparse pore/void architecture found in
% hydrogels. Voids have controllable size, size distribution, density, and
% shape roughness -- they are blob-like imperfect circles, NOT sharp cracks
% or perfect discs.
%
% -------------------------------------------------------------------------
% USAGE
%   [AtomsOut, BondsOut] = NetworkAddHeterogeneities(Atoms, Bonds, Domain, options)
%
% INPUTS
%   Atoms   [N x (5+MaxNbr)]  Atom array from network generator:
%               col 1    : atom ID  (consecutive integers)
%               col 2    : X position
%               col 3    : Y position
%               col 4    : Z position  (0 for 2-D networks)
%               col 5    : degree (number of bonds)
%               col 6..  : neighbor atom IDs (up to MaxNbr slots)
%
%   Bonds   [M x 4]           Bond array:
%               col 1    : bond ID
%               col 2    : atom ID 1
%               col 3    : atom ID 2
%               col 4    : bond length L
%
%   Domain  struct            Simulation domain (from NetworkGenSetup):
%               .xlo .xhi    x extents
%               .ylo .yhi    y extents
%               .min_node_sep  minimum node separation (used for defaults)
%
%   options struct            All fields are optional; defaults shown below.
%
% -------------------------------------------------------------------------
% OPTIONS
%
%  Void count / density
%   .density_mode    'count'          Total void count  (default)
%                    'area_fraction'  Drive coverage toward a target fraction
%   .n_voids         10               Number of voids  [density_mode='count']
%   .void_area_frac  0.10             Target areal coverage  [density_mode='area_fraction']
%
%  Void size
%   .radius_mean     3 * min_node_sep Mean void radius
%   .size_dist       'fixed'          All voids have radius_mean  (default)
%                    'gaussian'       Gaussian around radius_mean
%                    'exponential'    Exponential with mean = radius_mean
%   .radius_std      0.3*radius_mean  Std dev  [size_dist='gaussian']
%   .radius_min      0.5*radius_mean  Hard lower clip for sampled radii
%   .radius_max      3.0*radius_mean  Hard upper clip for sampled radii
%
%  Void shape  (rough imperfect circles via polar Fourier perturbation)
%   .shape_roughness 0.3              Fractional amplitude of boundary
%                                     perturbation. 0 = perfect circle,
%                                     ~0.5 = noticeably irregular blob,
%                                     ~0.8 = very ragged boundary.
%   .shape_n_modes   6                Number of angular Fourier modes.
%                                     Higher -> more complex, crinkled edge.
%                                     Recommended range: 3-12.
%
%  Void overlap
%   .void_overlap    true             true  : voids may overlap freely.
%                                     false : voids are kept separated by at
%                                             least bridge_width, guaranteeing
%                                             a strip of intact bonds between
%                                             every pair of voids.
%   .bridge_width    min_node_sep     Minimum edge-to-edge clearance enforced
%                                     when void_overlap = false.  Setting this
%                                     to a few multiples of min_node_sep
%                                     ensures at least one bond-width gap
%                                     between neighbouring voids.
%
%  Void center spatial distribution
%   .center_distribution  'random'   How void centers are distributed in space.
%
%     'random'    Uniform random placement -- centers drawn independently
%                 from a uniform distribution over the domain.  Classic
%                 Poisson-process-like behaviour; voids can cluster by chance.
%
%     'uniform'   Jittered grid -- the domain is divided into a grid of
%                 roughly equal cells (one per void) and each center is
%                 placed randomly within its own cell.  This prevents
%                 accidental clustering and gives even spatial coverage,
%                 like voids on a noisy lattice.
%
%     'clustered' Thomas cluster process -- a small number of invisible
%                 "parent" points are scattered first, then each void center
%                 is placed with a Gaussian offset from a randomly chosen
%                 parent.  Produces groups of nearby voids separated by
%                 larger void-free regions; good for modelling phase-
%                 separated or swelling-induced heterogeneity.
%
%   .n_cluster_parents  max(2, n_voids/4)
%                                     Number of cluster parent points
%                                     [center_distribution = 'clustered']
%   .cluster_spread     0.15*min(Lx,Ly)
%                                     Gaussian spread (std dev) of void
%                                     centers around each parent point.
%                                     Larger -> looser clusters.
%                                     [center_distribution = 'clustered']
%
%  Domain margin
%   .margin_frac     0.0              Fraction of radius_mean kept as margin
%                                     from domain boundary (0 = no margin)
%
%  Clamp regions  (defect-free bands at top and bottom of domain)
%   .clamp_frac      0.0              Fraction of total domain height (Ly) to
%                                     reserve as a clamp band at EACH end.
%                                     E.g. 0.1 protects the bottom 10% and top
%                                     10% of the domain from any void removal.
%                                     Void centers are never placed in these
%                                     bands, and atoms inside them are never
%                                     removed even if a void overlaps the edge.
%
% -------------------------------------------------------------------------
% OUTPUTS
%   AtomsOut  Pruned, renumbered Atoms array (same column layout as input)
%   BondsOut  Pruned, renumbered Bonds array (same column layout as input)
%
%   Atom and bond IDs are renumbered consecutively from 1.
%   Neighbor lists in AtomsOut are rebuilt to match BondsOut.
%
% -------------------------------------------------------------------------
% EXAMPLE  -- insert ~15% void coverage with irregular blobs
%
%   opts.density_mode    = 'area_fraction';
%   opts.void_area_frac  = 0.15;
%   opts.radius_mean     = 4 * Domain.min_node_sep;
%   opts.size_dist       = 'exponential';
%   opts.shape_roughness = 0.4;
%   opts.shape_n_modes   = 8;
%
%   [Atoms, Bonds] = NetworkAddHeterogeneities(Atoms, Bonds, Domain, opts);
%
% -------------------------------------------------------------------------
% NOTES
%  - R2016a compatible (no implicit expansion, no Statistics Toolbox).
%  - Network is assumed to be 2-D (Z column ignored for void checks).
%  - After removing atoms the Nvec array (if used) must be re-indexed:
%      keepBonds = ismember(old_BondIDs, BondsOut_original_IDs);
%      Nvec = Nvec(keepBonds);
%    where old_BondIDs are BondsOut(:,1) captured BEFORE calling this
%    function.  See the driver script for details.
% -------------------------------------------------------------------------

    %% ------------------------------------------------------------------ %%
    %                        0.  Parse options                              %
    %% ------------------------------------------------------------------ %%
    if nargin < 4,  options = struct();  end

    % Domain extents
    xlo = Domain.xlo;  xhi = Domain.xhi;
    ylo = Domain.ylo;  yhi = Domain.yhi;
    Lx  = xhi - xlo;
    Ly  = yhi - ylo;
    domain_area = Lx * Ly;

    % Default radius: 3x the minimum node separation is a natural scale
    if isfield(Domain, 'min_node_sep') && Domain.min_node_sep > 0
        default_r = 3.0 * Domain.min_node_sep;
    else
        default_r = 0.05 * min(Lx, Ly);
    end

    options = set_default(options, 'density_mode',    'count');
    options = set_default(options, 'n_voids',         10);
    options = set_default(options, 'void_area_frac',  0.10);
    options = set_default(options, 'radius_mean',     default_r);
    options = set_default(options, 'size_dist',       'fixed');
    % Defer std/min/max defaults until after radius_mean is confirmed present
    options = set_default(options, 'radius_std',      0.30 * options.radius_mean);
    options = set_default(options, 'radius_min',      0.50 * options.radius_mean);
    options = set_default(options, 'radius_max',      3.00 * options.radius_mean);
    options = set_default(options, 'shape_roughness',     0.30);
    options = set_default(options, 'shape_n_modes',       6);
    options = set_default(options, 'void_overlap',        true);
    options = set_default(options, 'bridge_width',        Domain.min_node_sep);
    options = set_default(options, 'center_distribution', 'random');
    options = set_default(options, 'margin_frac',         0.0);
    options = set_default(options, 'clamp_frac',          0.0);
    % n_cluster_parents and cluster_spread depend on n_voids / domain size,
    % so their defaults are set later, after n_voids is determined.

    r_mean    = options.radius_mean;
    roughness = options.shape_roughness;
    n_modes   = max(1, round(options.shape_n_modes));

    %% ------------------------------------------------------------------ %%
    %              1.  Determine number of voids to place                   %
    %% ------------------------------------------------------------------ %%
    if strcmpi(options.density_mode, 'area_fraction')
        mean_void_area = pi * r_mean^2;
        n_voids = max(1, round(options.void_area_frac * domain_area / mean_void_area));
        fprintf('   [Heterogeneities] area_fraction mode: targeting %.1f%% coverage -> %d voids\n', ...
            options.void_area_frac * 100, n_voids);
    else
        n_voids = max(1, round(options.n_voids));
        fprintf('   [Heterogeneities] count mode: %d voids requested\n', n_voids);
    end

    %% ------------------------------------------------------------------ %%
    %                  2.  Sample void radii                                %
    %% ------------------------------------------------------------------ %%
    radii = sample_radii(n_voids, options);   % [n_voids x 1]

    %% ------------------------------------------------------------------ %%
    %                  3.  Place void centers                               %
    %% ------------------------------------------------------------------ %%
    margin = options.margin_frac * r_mean;

    x_lo = xlo + margin;   x_hi = xhi - margin;
    y_lo = ylo + margin;   y_hi = yhi - margin;

    % Safety: if margin made the range invalid, revert to full domain
    if x_hi <= x_lo,  x_lo = xlo;  x_hi = xhi;  end
    if y_hi <= y_lo,  y_lo = ylo;  y_hi = yhi;  end

    eff_Lx = x_hi - x_lo;
    eff_Ly = y_hi - y_lo;

    % --- Clamp bands: bottom and top strips that must remain defect-free ---
    clamp_h = options.clamp_frac * Ly;   % absolute height of each clamp band
    clamp_ylo = ylo + clamp_h;           % lower edge of the active (non-clamped) zone
    clamp_yhi = yhi - clamp_h;           % upper edge of the active (non-clamped) zone
    if clamp_yhi <= clamp_ylo
        warning('NetworkAddHeterogeneities: clamp_frac=%.2g leaves no active zone; ignoring clamps.', ...
            options.clamp_frac);
        clamp_h   = 0;
        clamp_ylo = ylo;
        clamp_yhi = yhi;
    end

    % Restrict void center placement to the active zone (tighter of margin and clamp)
    y_lo = max(y_lo, clamp_ylo);
    y_hi = min(y_hi, clamp_yhi);
    if y_hi <= y_lo,  y_lo = clamp_ylo;  y_hi = clamp_yhi;  end
    eff_Ly = y_hi - y_lo;   % recompute after clamping

    % --- Set cluster defaults now that n_voids and domain size are known ---
    options = set_default(options, 'n_cluster_parents', max(2, round(n_voids / 4)));
    options = set_default(options, 'cluster_spread',    0.15 * min(eff_Lx, eff_Ly));

    % -----------------------------------------------------------------------
    %  3a.  Generate candidate centers according to center_distribution
    % -----------------------------------------------------------------------
    switch lower(options.center_distribution)

        case 'random'
            % Uniform independent placement -- classic Poisson-like behaviour.
            % Voids can cluster by chance; no spatial regularity enforced.
            cand_x = x_lo + eff_Lx * rand(n_voids, 1);
            cand_y = y_lo + eff_Ly * rand(n_voids, 1);

        case 'uniform'
            % Jittered grid -- domain divided into n_voids roughly-square cells;
            % one center placed randomly inside each cell.  Prevents accidental
            % clustering and gives even spatial coverage.
            nx_g = max(1, round(sqrt(n_voids * eff_Lx / eff_Ly)));
            ny_g = max(1, ceil(n_voids / nx_g));
            cell_w = eff_Lx / nx_g;
            cell_h = eff_Ly / ny_g;

            % Build flat list of (col, row) cell indices and shuffle
            [col_idx, row_idx] = meshgrid(0:nx_g-1, 0:ny_g-1);
            col_flat = col_idx(:);
            row_flat = row_idx(:);
            n_cells  = numel(col_flat);
            perm     = randperm(n_cells);
            col_flat = col_flat(perm);
            row_flat = row_flat(perm);

            % Take the first n_voids cells; place center uniformly within each
            use = min(n_voids, n_cells);
            cand_x = zeros(n_voids, 1);
            cand_y = zeros(n_voids, 1);
            for v = 1:use
                cand_x(v) = x_lo + (col_flat(v) + rand) * cell_w;
                cand_y(v) = y_lo + (row_flat(v) + rand) * cell_h;
            end
            if use < n_voids
                % More voids than cells: fill remainder with random placement
                cand_x(use+1:end) = x_lo + eff_Lx * rand(n_voids - use, 1);
                cand_y(use+1:end) = y_lo + eff_Ly * rand(n_voids - use, 1);
            end

        case 'clustered'
            % Thomas cluster process -- scatter a few parent points, then place
            % each void center as a Gaussian offset from a random parent.
            % Produces groups of nearby voids separated by larger empty regions.
            n_par  = options.n_cluster_parents;
            spread = options.cluster_spread;

            % Place parent points uniformly in domain
            par_x = x_lo + eff_Lx * rand(n_par, 1);
            par_y = y_lo + eff_Ly * rand(n_par, 1);

            % Assign each void to a random parent and offset by Gaussian
            cand_x = zeros(n_voids, 1);
            cand_y = zeros(n_voids, 1);
            for v = 1:n_voids
                p      = randi(n_par);
                % Clamp to domain with a while loop (rare for reasonable spread)
                for attempt = 1:50
                    cx_try = par_x(p) + spread * randn;
                    cy_try = par_y(p) + spread * randn;
                    if cx_try >= x_lo && cx_try <= x_hi && ...
                       cy_try >= y_lo && cy_try <= y_hi
                        break;
                    end
                end
                % If still out of bounds after retries, clamp hard
                cand_x(v) = max(x_lo, min(x_hi, cx_try));
                cand_y(v) = max(y_lo, min(y_hi, cy_try));
            end

        otherwise
            warning(['NetworkAddHeterogeneities: unknown center_distribution ' ...
                     '"%s"; falling back to random.'], options.center_distribution);
            cand_x = x_lo + eff_Lx * rand(n_voids, 1);
            cand_y = y_lo + eff_Ly * rand(n_voids, 1);
    end

    fprintf('   [Heterogeneities] center_distribution = %s\n', options.center_distribution);

    % -----------------------------------------------------------------------
    %  3b.  Enforce overlap / bridge constraint
    %
    %  When void_overlap = false, two voids must not overlap AND must be
    %  separated by at least bridge_width at their closest edges, ensuring
    %  an intact strip of bonds runs between every pair.
    %
    %  IMPORTANT: The actual void boundary extends up to (1+roughness)*R
    %  in the worst case (Fourier perturbation).  We use this worst-case
    %  effective radius in the separation test so that rough boundaries
    %  cannot sneak through the gap check.
    %
    %  Effective edge-to-edge gap = d_centers - eff_r_i - eff_r_j >= bridge_width
    %  => required: d_centers >= eff_r_i + eff_r_j + bridge_width
    % -----------------------------------------------------------------------
    centers    = zeros(n_voids, 2);
    eff_radii  = radii * (1 + roughness);   % worst-case boundary extent per void
    placed     = 0;

    % Normalise void_overlap: accept logical, numeric, or string 'true'/'false'
    % Guards against the common mistake of passing the string 'false'
    % instead of logical false -- a non-empty string is always truthy in MATLAB.
    vo = options.void_overlap;
    if ischar(vo) || isstring(vo)
        vo = ~(strcmpi(vo,'false') || strcmp(vo,'0'));
    end
    vo = logical(vo);

    if vo
        % --- Overlap allowed: accept all candidates as-is ---
        centers = [cand_x, cand_y];
        placed  = n_voids;

    else
        % --- Overlap forbidden: test each candidate against already-placed voids
        bridge           = max(0, options.bridge_width);
        max_tries_per_void = 200;   % independent per-void retry budget

        for v = 1:n_voids
            accepted = false;

            % Try the pre-generated candidate first; then re-sample randomly
            for attempt = 1 : max_tries_per_void
                if attempt == 1
                    cx_try = cand_x(v);
                    cy_try = cand_y(v);
                else
                    cx_try = x_lo + eff_Lx * rand;
                    cy_try = y_lo + eff_Ly * rand;
                end

                ok = true;
                for prev = 1:placed
                    % Use effective (roughness-inflated) radii so that the
                    % rough boundaries of both voids are guaranteed not to meet
                    d_edge = sqrt((cx_try - centers(prev,1))^2 + ...
                                  (cy_try - centers(prev,2))^2) ...
                             - eff_radii(v) - eff_radii(prev);
                    if d_edge < bridge
                        ok = false;
                        break;
                    end
                end

                if ok
                    placed = placed + 1;
                    centers(placed, 1) = cx_try;
                    centers(placed, 2) = cy_try;
                    accepted = true;
                    break;
                end
            end

            if ~accepted
                warning(['NetworkAddHeterogeneities: could not place void %d ' ...
                         'with bridge_width=%.2g after %d attempts; ' ...
                         'stopping at %d placed voids. ' ...
                         'Consider reducing n_voids, radius_mean, shape_roughness, ' ...
                         'or bridge_width.'], ...
                         v, bridge, max_tries_per_void, placed);
                break;
            end
        end

        n_voids   = placed;
        radii     = radii(1:placed);
        eff_radii = eff_radii(1:placed);
        centers   = centers(1:placed, :);
    end

    fprintf('   [Heterogeneities] Placed %d voids  (void_overlap=%d, bridge_width=%.2g)\n', ...
        n_voids, vo, options.bridge_width);

    %% ------------------------------------------------------------------ %%
    %          4.  Generate per-void boundary roughness parameters          %
    %                                                                       %
    %  The void boundary in polar coordinates is:                          %
    %    r(theta) = R * [1 + A * sum_k( w_k * cos(k*theta + phi_k) )]     %
    %                                                                       %
    %  where A = shape_roughness, w_k are amplitude weights (normalized,   %
    %  decaying as 1/k for low-frequency blob-like shapes), and phi_k are  %
    %  independent uniform random phases.                                  %
    %                                                                       %
    %  Mode indices run from 2..n_modes+1 (k=1 is skipped: it would just   %
    %  shift the apparent center rather than deform the shape).            %
    %% ------------------------------------------------------------------ %%
    mode_weights = zeros(n_voids, n_modes);
    mode_phases  = zeros(n_voids, n_modes);

    for v = 1:n_voids
        % Base amplitudes decay as 1/k (preference for large-scale deformation)
        k_vals   = 2 : (n_modes + 1);              % skip k=1
        raw_amps = (1 ./ k_vals) .* (0.5 + rand(1, n_modes));
        % Normalize so weights sum to 1 (roughness parameter controls total amplitude)
        mode_weights(v, :) = raw_amps / sum(raw_amps);
        mode_phases(v, :)  = 2 * pi * rand(1, n_modes);
    end

    %% ------------------------------------------------------------------ %%
    %              5.  Classify atoms as inside or outside voids            %
    %% ------------------------------------------------------------------ %%
    natom   = size(Atoms, 1);
    in_void = false(natom, 1);

    atom_x = Atoms(:, 2);
    atom_y = Atoms(:, 3);

    for v = 1:n_voids
        cx_v    = centers(v, 1);
        cy_v    = centers(v, 2);
        R_v     = radii(v);
        wts     = mode_weights(v, :);   % [1 x n_modes]
        phis    = mode_phases(v, :);    % [1 x n_modes]
        k_vals  = 2 : (n_modes + 1);   % [1 x n_modes]

        dx = atom_x - cx_v;
        dy = atom_y - cy_v;
        dist2_all = dx .* dx + dy .* dy;

        % Bounding-circle pre-filter: only examine atoms within the
        % maximum possible void radius to avoid expensive computation
        outer_R   = (1.0 + roughness) * R_v;
        candidate = find(dist2_all < outer_R * outer_R & ~in_void);

        if isempty(candidate),  continue;  end

        % For candidate atoms compute local void radius at their angle
        theta_c = atan2(dy(candidate), dx(candidate));  % [nc x 1]

        % Vectorised Fourier sum over modes
        % perturbation(i) = sum_m( wts(m) * cos(k_vals(m)*theta_c(i) + phis(m)) )
        % Build angle-mode matrix: [nc x n_modes]
        k_mat   = repmat(k_vals,  numel(candidate), 1);   % [nc x nm]
        phi_mat = repmat(phis,    numel(candidate), 1);   % [nc x nm]
        w_mat   = repmat(wts,     numel(candidate), 1);   % [nc x nm]
        th_mat  = repmat(theta_c, 1, n_modes);            % [nc x nm]

        perturbation = sum(w_mat .* cos(k_mat .* th_mat + phi_mat), 2);  % [nc x 1]

        r_local = R_v * max(0.1, 1 + roughness * perturbation);  % [nc x 1]
        % (max ensures r_local stays positive even for strong roughness)

        dist_c = sqrt(dist2_all(candidate));   % [nc x 1]

        in_void(candidate(dist_c < r_local)) = true;
    end

    % Protect atoms inside the clamp bands -- they must never be removed
    if clamp_h > 0
        in_clamp         = (atom_y <= clamp_ylo) | (atom_y >= clamp_yhi);
        in_void(in_clamp) = false;
    end

    n_removed = sum(in_void);
    fprintf('   [Heterogeneities] Marking %d / %d atoms for removal (%.1f%%)\n', ...
        n_removed, natom, 100.0 * n_removed / max(natom, 1));

    %% ------------------------------------------------------------------ %%
    %                  6.  Remove atoms inside voids                        %
    %% ------------------------------------------------------------------ %%
    AtomsOut = Atoms(~in_void, :);

    if isempty(AtomsOut)
        warning(['NetworkAddHeterogeneities: ALL atoms were removed. ' ...
                 'Void radii or count may be too large for this domain.']);
        AtomsOut = zeros(0, size(Atoms, 2));
        BondsOut = zeros(0, 4);
        return;
    end

    surviving_ids = AtomsOut(:, 1);   % IDs still present

    %% ------------------------------------------------------------------ %%
    %               7.  Remove bonds incident to removed atoms              %
    %% ------------------------------------------------------------------ %%
    if isempty(Bonds)
        BondsOut = zeros(0, 4);
    else
        % Keep a bond only if BOTH endpoints survived
        keep1    = ismember(Bonds(:, 2), surviving_ids);
        keep2    = ismember(Bonds(:, 3), surviving_ids);
        BondsOut = Bonds(keep1 & keep2, :);
    end

    n_bonds_rm = size(Bonds, 1) - size(BondsOut, 1);
    fprintf('   [Heterogeneities] Removed %d / %d bonds\n', ...
        n_bonds_rm, size(Bonds, 1));

    %% ------------------------------------------------------------------ %%
    %            8.  Renumber atom IDs and bond endpoints                   %
    %% ------------------------------------------------------------------ %%
    natom_new     = size(AtomsOut, 1);
    old_ids       = AtomsOut(:, 1);          % original (possibly non-consecutive) IDs
    new_ids       = (1 : natom_new)';        % new consecutive IDs
    AtomsOut(:,1) = new_ids;

    if ~isempty(BondsOut)
        % Remap bond endpoints from old IDs to new IDs
        [~, loc1]    = ismember(BondsOut(:, 2), old_ids);
        [~, loc2]    = ismember(BondsOut(:, 3), old_ids);
        BondsOut(:,2) = new_ids(loc1);
        BondsOut(:,3) = new_ids(loc2);
        BondsOut(:,1) = (1 : size(BondsOut, 1))';   % renumber bond IDs
    end

    %% ------------------------------------------------------------------ %%
    %               9.  Rebuild neighbor lists in AtomsOut                  %
    %% ------------------------------------------------------------------ %%
    Max_peratom_bond = size(AtomsOut, 2) - 5;
    if Max_peratom_bond < 1,  Max_peratom_bond = 1;  end

    AtomsOut(:, 5)                        = 0;   % reset degree
    AtomsOut(:, 6 : 5 + Max_peratom_bond) = 0;   % clear neighbor slots

    for k = 1 : size(BondsOut, 1)
        ii  = BondsOut(k, 2);
        jj  = BondsOut(k, 3);

        nb_i = AtomsOut(ii, 5) + 1;
        if nb_i <= Max_peratom_bond
            AtomsOut(ii, 5)        = nb_i;
            AtomsOut(ii, 5 + nb_i) = jj;
        end

        nb_j = AtomsOut(jj, 5) + 1;
        if nb_j <= Max_peratom_bond
            AtomsOut(jj, 5)        = nb_j;
            AtomsOut(jj, 5 + nb_j) = ii;
        end
    end

    %% ------------------------------------------------------------------ %%
    %                        10.  Summary                                   %
    %% ------------------------------------------------------------------ %%
    fprintf(['   [Heterogeneities] Complete. ' ...
             '%d atoms remain (-%d), %d bonds remain (-%d).\n'], ...
        natom_new,         n_removed, ...
        size(BondsOut, 1), n_bonds_rm);

end   % end main function


%% ======================================================================= %%
%                           LOCAL HELPER FUNCTIONS                          %
%% ======================================================================= %%

function radii = sample_radii(n, opts)
% SAMPLE_RADII  Draw n void radii from the requested distribution.
%
% Supported distributions (opts.size_dist):
%   'fixed'       All radii equal opts.radius_mean
%   'gaussian'    Normal(radius_mean, radius_std), clipped to [radius_min, radius_max]
%   'exponential' Exponential with mean=radius_mean, shifted to start at radius_min,
%                 clipped at radius_max.  No Statistics Toolbox required.

    r_mean = opts.radius_mean;
    r_min  = opts.radius_min;
    r_max  = opts.radius_max;

    switch lower(opts.size_dist)

        case 'fixed'
            radii = r_mean * ones(n, 1);

        case 'gaussian'
            r_std = opts.radius_std;
            radii = r_mean + r_std * randn(n, 1);
            radii = max(r_min, min(r_max, radii));

        case 'exponential'
            % Inverse-transform sampling (no toolbox needed):
            %   If U ~ Uniform(0,1), then X = -mu * log(U) ~ Exp(mu)
            % We shift by r_min so the distribution starts there.
            adjusted_mean = max(r_mean - r_min, 1e-10);
            U     = rand(n, 1);
            U     = max(U, 1e-15);          % avoid log(0)
            radii = r_min + (-adjusted_mean * log(U));
            radii = min(radii, r_max);      % cap the tail

        otherwise
            warning('NetworkAddHeterogeneities:sample_radii:unknownDist', ...
                'Unknown size_dist "%s"; defaulting to fixed.', opts.size_dist);
            radii = r_mean * ones(n, 1);
    end
end


function s = set_default(s, field, val)
% SET_DEFAULT  Assign val to s.field only if the field does not already exist.
    if ~isfield(s, field)
        s.(field) = val;
    end
end