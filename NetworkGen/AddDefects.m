function [Atoms, Bonds, Nvec] = AddDefects(obj, Atoms, Bonds, Nvec)
% -------------------------------------------------------------------------
% AddDefects
%   Insert heterogeneous void defects into a polymer network.
%
%   Reads all settings from obj.defect, obj.domain, and obj.flags.
%   Void shape, placement, size distribution, and sparse/foam mode are
%   all controlled by obj.defect fields (see network.m for defaults).
%
%   Three sequential passes:
%     Pass 1 — Primary void placement  (removes atoms inside voids)
%     Pass 2 — Density thinning        (thins over-dense bridges/patches)
%     Pass 3 — Constriction bridging   (adds single bonds across narrow gaps)
%
%   NOTE: This function does NOT prune isolated atoms, remove disconnected
%   clusters, or renumber IDs. All cleanup is deferred to CleanupNetwork.
%
% INPUT
%   obj   : network object (reads obj.defect, obj.domain, obj.flags)
%   Atoms : [N x (5+MaxNbr)]  atom array
%   Bonds : [M x 5]           bond array  [id | i | j | L0 | type]
%   Nvec  : [M x ...]         per-bond quantity array (e.g. Kuhn-N values)
%
% OUTPUT
%   Atoms : atom array after all passes
%   Bonds : bond array after all passes (not yet renumbered)
%   Nvec  : Nvec filtered/extended to match surviving/added bonds
% -------------------------------------------------------------------------

    %% ------------------------------------------------------------------
    %  0.  Early exit if defects are disabled
    %% ------------------------------------------------------------------
    if ~obj.flags.idefect
        obj.log.print('   [AddDefects] Skipped (flags.idefect = false)\n');
        return;
    end

    %% ------------------------------------------------------------------
    %  1.  Unpack settings from obj
    %% ------------------------------------------------------------------
    d = obj.defect;

    xlo = obj.domain.xlo;  xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;  yhi = obj.domain.yhi;
    Lx  = xhi - xlo;
    Ly  = yhi - ylo;
    domain_area = Lx * Ly;

    r_mean    = d.radius_mean;
    r_std     = d.radius_std;
    r_min     = d.radius_min;
    r_max     = d.radius_max;
    roughness = d.shape_roughness;
    n_modes   = max(1, round(d.shape_n_modes));
    wall_t    = d.wall_thickness;
    do_sparse = logical(d.sparse_network);

    %% ------------------------------------------------------------------
    %  2.  Determine number of voids
    %% ------------------------------------------------------------------
    if strcmpi(d.density_mode, 'area_frac')
        mean_void_area = pi * r_mean^2;
        n_voids = max(1, round(d.void_area_frac * domain_area / mean_void_area));
        obj.log.print('   [AddDefects] area_frac mode: targeting %.1f%% coverage -> %d voids\n', ...
            d.void_area_frac * 100, n_voids);
    else
        n_voids = max(1, round(d.n_voids));
        obj.log.print('   [AddDefects] count mode: %d voids requested\n', n_voids);
    end

    %% ------------------------------------------------------------------
    %  3.  Sample void radii
    %% ------------------------------------------------------------------
    radii = local_sample_radii(n_voids, d.size_dist, r_mean, r_std, r_min, r_max);

    %% ------------------------------------------------------------------
    %  4.  Compute placement zone (margin + clamp bands)
    %% ------------------------------------------------------------------
    margin  = d.margin_frac * r_mean;
    x_lo    = xlo + margin;   x_hi = xhi - margin;
    y_lo    = ylo + margin;   y_hi = yhi - margin;

    % Safety: if margin collapsed the range, revert to full domain
    if x_hi <= x_lo,  x_lo = xlo;  x_hi = xhi;  end
    if y_hi <= y_lo,  y_lo = ylo;  y_hi = yhi;  end

    % Clamp bands (defect-free strips at top and bottom)
    clamp_h   = d.clamp_thickness * Ly;
    clamp_ylo = ylo + clamp_h;
    clamp_yhi = yhi - clamp_h;

    if clamp_yhi <= clamp_ylo
        warning('AddDefects: clamp_thickness=%.2g leaves no active zone; clamps ignored.', ...
            d.clamp_thickness);
        clamp_h   = 0;
        clamp_ylo = ylo;
        clamp_yhi = yhi;
    end

    % Restrict void placement to the tighter of margin and clamp
    y_lo = max(y_lo, clamp_ylo);
    y_hi = min(y_hi, clamp_yhi);
    if y_hi <= y_lo
        y_lo = clamp_ylo;
        y_hi = clamp_yhi;
    end

    eff_Lx = x_hi - x_lo;
    eff_Ly = y_hi - y_lo;

    %% ------------------------------------------------------------------
    %  5.  Place void centers
    %% ------------------------------------------------------------------
    switch lower(d.center_distribution)

        case 'random'
            cand_x = x_lo + eff_Lx * rand(n_voids, 1);
            cand_y = y_lo + eff_Ly * rand(n_voids, 1);

        case 'uniform'
            nx_g   = max(1, round(sqrt(n_voids * eff_Lx / eff_Ly)));
            ny_g   = max(1, ceil(n_voids / nx_g));
            cell_w = eff_Lx / nx_g;
            cell_h = eff_Ly / ny_g;
            [gx, gy] = meshgrid( ...
                x_lo + cell_w * (0.5 + (0:nx_g-1)), ...
                y_lo + cell_h * (0.5 + (0:ny_g-1)));
            pts    = [gx(:), gy(:)];
            pts    = pts(randperm(size(pts,1)), :);
            pts    = pts(1:min(n_voids, end), :);
            cand_x = pts(:,1) + cell_w * 0.3 * (rand(size(pts,1),1) - 0.5);
            cand_y = pts(:,2) + cell_h * 0.3 * (rand(size(pts,1),1) - 0.5);
            cand_x = max(x_lo, min(x_hi, cand_x));
            cand_y = max(y_lo, min(y_hi, cand_y));

        case 'clustered'
            np     = max(1, round(d.n_cluster_parents));
            spread = max(0, d.cluster_spread);
            px     = x_lo + eff_Lx * rand(np, 1);
            py     = y_lo + eff_Ly * rand(np, 1);
            assign = randi(np, n_voids, 1);
            cand_x = px(assign) + spread * randn(n_voids, 1);
            cand_y = py(assign) + spread * randn(n_voids, 1);
            cand_x = max(x_lo, min(x_hi, cand_x));
            cand_y = max(y_lo, min(y_hi, cand_y));

        case 'manual'
            % Use manually-specified void center coordinates
            if ~strcmpi(d.density_mode, 'count')
                error('AddDefects: manual center_distribution only works with density_mode=''count''');
            end
            if isempty(d.manual_centers)
                error('AddDefects: manual center_distribution requires manual_centers to be set');
            end
            if size(d.manual_centers, 1) ~= n_voids || size(d.manual_centers, 2) ~= 2
                error('AddDefects: manual_centers must be [%d x 2] array but got [%d x %d]', ...
                    n_voids, size(d.manual_centers, 1), size(d.manual_centers, 2));
            end
            cand_x = d.manual_centers(:, 1);
            cand_y = d.manual_centers(:, 2);
            % Clamp to domain bounds
            cand_x = max(xlo, min(xhi, cand_x));
            cand_y = max(ylo, min(yhi, cand_y));
            obj.log.print('   [AddDefects] Using manual coordinates for %d voids\n', n_voids);

        otherwise
            warning('AddDefects: unknown center_distribution "%s"; using random.', ...
                d.center_distribution);
            cand_x = x_lo + eff_Lx * rand(n_voids, 1);
            cand_y = y_lo + eff_Ly * rand(n_voids, 1);
    end

    %% ------------------------------------------------------------------
    %  6.  Resolve overlaps (bridge_width enforcement)
    %% ------------------------------------------------------------------
    eff_radii = radii;
    centers   = zeros(n_voids, 2);
    placed    = 0;
    vo        = logical(d.void_overlap);

    if vo
        centers   = [cand_x, cand_y];
        placed    = n_voids;
    else
        bridge           = max(0, d.bridge_width);
        max_tries_per_void = 200;

        for v = 1:n_voids
            accepted = false;
            for attempt = 1:max_tries_per_void
                if attempt == 1
                    cx_try = cand_x(v);
                    cy_try = cand_y(v);
                else
                    cx_try = x_lo + eff_Lx * rand;
                    cy_try = y_lo + eff_Ly * rand;
                end

                ok = true;
                for prev = 1:placed
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
                    centers(placed,1) = cx_try;
                    centers(placed,2) = cy_try;
                    accepted = true;
                    break;
                end
            end

            if ~accepted
                warning(['AddDefects: could not place void %d with ' ...
                         'bridge_width=%.2g after %d attempts; ' ...
                         'stopping at %d placed voids.'], ...
                         v, bridge, max_tries_per_void, placed);
                break;
            end
        end

        n_voids   = placed;
        radii     = radii(1:placed);
        eff_radii = eff_radii(1:placed);
        centers   = centers(1:placed, :);
    end

    obj.log.print('   [AddDefects] Placed %d voids  (void_overlap=%d)\n', n_voids, vo);

    %% ------------------------------------------------------------------
    %  7.  Generate per-void boundary roughness
    %% ------------------------------------------------------------------
    mode_weights = zeros(n_voids, n_modes);
    mode_phases  = zeros(n_voids, n_modes);
    k_vals       = 2:(n_modes + 1);

    for v = 1:n_voids
        raw_amps           = (1 ./ k_vals) .* (0.5 + rand(1, n_modes));
        mode_weights(v,:)  = raw_amps / sum(raw_amps);
        mode_phases(v,:)   = 2 * pi * rand(1, n_modes);
    end

    %% ------------------------------------------------------------------
    %  8.  Classify atoms as inside void / near boundary
    %% ------------------------------------------------------------------
    natom        = size(Atoms, 1);
    in_void       = false(natom, 1);
    near_boundary = false(natom, 1);

    % Atoms column layout (post-refactor): [ID | molID | X | Y | Z | deg | nbrs]
    atom_x = Atoms(:, 3);
    atom_y = Atoms(:, 4);

    for v = 1:n_voids
        cx_v = centers(v, 1);
        cy_v = centers(v, 2);
        R_v  = radii(v);
        wts  = mode_weights(v, :);
        phis = mode_phases(v, :);

        dx       = atom_x - cx_v;
        dy       = atom_y - cy_v;
        dist2    = dx .* dx + dy .* dy;

        outer_R   = (1.0 + roughness) * R_v + wall_t;
        candidate = find(dist2 < outer_R * outer_R);
        if isempty(candidate),  continue;  end

        theta_c  = atan2(dy(candidate), dx(candidate));

        k_mat    = repmat(k_vals,     numel(candidate), 1);
        phi_mat  = repmat(phis,       numel(candidate), 1);
        w_mat    = repmat(wts,        numel(candidate), 1);
        th_mat   = repmat(theta_c,    1, n_modes);

        perturb  = sum(w_mat .* cos(k_mat .* th_mat + phi_mat), 2);
        r_local  = R_v * max(0.1, 1 + roughness * perturb);

        dist_c   = sqrt(dist2(candidate));
        inside   = dist_c < r_local;

        in_void(candidate(inside)) = true;

        if do_sparse
            near = ~inside & (dist_c < r_local + wall_t);
            near_boundary(candidate(near)) = true;
        end
    end

    % Protect clamp-band atoms (y-boundaries)
    if clamp_h > 0
        in_clamp               = (atom_y <= clamp_ylo) | (atom_y >= clamp_yhi);
        in_void(in_clamp)       = false;
        near_boundary(in_clamp) = true;
    end

    % Protect x-boundary atoms in sparse mode.
    % Only applied for 'fixed' boundary — periodic networks wrap around
    % and need no boundary protection on x-edges.
    if do_sparse && strcmpi(obj.domain.boundary, 'fixed')
        in_x_wall                = (atom_x <= xlo + wall_t) | (atom_x >= xhi - wall_t);
        in_void(in_x_wall)       = false;
        near_boundary(in_x_wall) = true;
    end

    % Final removal mask
    if do_sparse
        remove_mask = in_void | ~near_boundary;
        obj.log.print('   [AddDefects] sparse_network mode  (wall_thickness=%.2g)\n', wall_t);
    else
        remove_mask = in_void;
    end

    n_removed = sum(remove_mask);
    obj.log.print('   [AddDefects] Marking %d / %d atoms for removal (%.1f%%)\n', ...
        n_removed, natom, 100.0 * n_removed / max(natom, 1));

    %% ------------------------------------------------------------------
    %  9.  Remove atoms inside voids
    %% ------------------------------------------------------------------
    Atoms = Atoms(~remove_mask, :);

    if isempty(Atoms)
        warning(['AddDefects: ALL atoms were removed. ' ...
                 'Void radii or count may be too large for this domain.']);
        Bonds = zeros(0, size(Bonds, 2));
        Nvec  = [];
        return;
    end

    surviving_ids = Atoms(:, 1);

    %% ------------------------------------------------------------------
    %  10.  Remove bonds whose endpoints were removed
    %% ------------------------------------------------------------------
    if isempty(Bonds)
        Bonds = zeros(0, size(Bonds, 2));
        Nvec  = [];
        obj.log.print('   [AddDefects] No bonds to filter\n');
        return;
    end

    keep1    = ismember(Bonds(:, 2), surviving_ids);
    keep2    = ismember(Bonds(:, 3), surviving_ids);
    keep     = keep1 & keep2;

    n_bonds_rm = sum(~keep);
    Bonds = Bonds(keep, :);
    obj.log.print('   [AddDefects] Removed %d / %d bonds incident to void atoms\n', ...
        n_bonds_rm, n_bonds_rm + size(Bonds, 1));

    %% ------------------------------------------------------------------
    %  11.  Filter Nvec to match surviving bonds
    %% ------------------------------------------------------------------
    if ~isempty(Nvec)
        Nvec = Nvec(keep, :);
    end

    obj.log.print(['   [AddDefects] Pass 1 done. ' ...
             '%d atoms remain (-%d), %d bonds remain (-%d).\n'], ...
        size(Atoms, 1),  n_removed, ...
        size(Bonds, 1),  n_bonds_rm);

    %% ------------------------------------------------------------------
    %  12.  Pass 2 — Secondary density thinning
    %% ------------------------------------------------------------------
    if isfield(d, 'thinning') && logical(d.thinning) && ~isempty(Atoms) && ~isempty(Bonds)
        [Atoms, Bonds, Nvec] = local_thin_network(obj, Atoms, Bonds, Nvec);
    end

    %% ------------------------------------------------------------------
    %  13.  Pass 3 — Constriction bridging
    %        Adds single bonds across narrow void constrictions to prevent
    %        unrealistically long unbroken void channels.
    %% ------------------------------------------------------------------
    if isfield(d, 'bridging') && logical(d.bridging) && ~isempty(Atoms) && ~isempty(Bonds)
        [Atoms, Bonds, Nvec] = local_bridge_constrictions(obj, Atoms, Bonds, Nvec);
    end

end   % end main function


%% =======================================================================
%  LOCAL HELPER:  local_sample_radii
%% =======================================================================
function radii = local_sample_radii(n, size_dist, r_mean, r_std, r_min, r_max)
% Draw n void radii from the requested distribution.

    switch lower(size_dist)

        case 'fixed'
            radii = r_mean * ones(n, 1);

        case 'gaussian'
            radii = r_mean + r_std * randn(n, 1);
            radii = max(r_min, min(r_max, radii));

        case 'exponential'
            adjusted_mean = max(r_mean - r_min, 1e-10);
            U     = rand(n, 1);
            U     = max(U, 1e-15);
            radii = r_min + (-adjusted_mean * log(U));
            radii = min(radii, r_max);

        otherwise
            warning('AddDefects:sample_radii:unknownDist', ...
                'Unknown size_dist "%s"; defaulting to fixed.', size_dist);
            radii = r_mean * ones(n, 1);
    end
end


%% =======================================================================
%  LOCAL HELPER:  local_thin_network
%
%  Pass 2 — probabilistic atom removal in over-dense regions.
%
%  For each atom, estimates local density by counting neighbours within
%  thinning_radius. Atoms denser than average get a removal probability
%  proportional to how over-dense they are. Atoms at or below average
%  density are left alone.
%
%  Reads from obj.defect:
%    thinning_radius      - neighbourhood radius for density estimate
%                           0 = auto (2.5 × avg node spacing)
%    thinning_target_frac - target keep fraction in over-dense regions
%                           (default 0.4)
%    thinning_min_keep    - minimum keep fraction anywhere (default 0.1)
%% =======================================================================
function [Atoms, Bonds, Nvec] = local_thin_network(obj, Atoms, Bonds, Nvec)

    d = obj.defect;

    rho_atom    = obj.architecture.rho_atom;
    avg_spacing = 1.0 / sqrt(max(rho_atom, 1e-10));

    b = obj.domain.b;
    r_thin = 2.5 * avg_spacing;            % auto default in real units
    if isfield(d, 'thinning_radius') && d.thinning_radius > 0
        r_thin = d.thinning_radius * b;    % user value in units of b
    end

    target_frac = 0.4;
    if isfield(d, 'thinning_target_frac')
        target_frac = d.thinning_target_frac;
    end

    min_keep = 0.1;
    if isfield(d, 'thinning_min_keep')
        min_keep = d.thinning_min_keep;
    end

    obj.log.print('   [AddDefects:Thin] Pass 2 — density thinning\n');
    obj.log.print('   [AddDefects:Thin] r_thin=%.1f  target_frac=%.3f  min_keep=%.3f\n', ...
        r_thin, target_frac, min_keep);

    natom_before = size(Atoms, 1);
    atom_ids     = Atoms(:, 1);
    % Atoms column layout: [ID | molID | X | Y | Z | deg | nbrs]
    ax           = Atoms(:, 3);
    ay           = Atoms(:, 4);
    natom        = numel(atom_ids);

    % Expected atoms in a circle of radius r_thin under average density
    expected_in_r = max(rho_atom * pi * r_thin^2, 1.0);

    % Local density per atom
    local_density = zeros(natom, 1);
    for i = 1:natom
        dx2 = (ax - ax(i)).^2 + (ay - ay(i)).^2;
        local_density(i) = (sum(dx2 <= r_thin^2) - 1) / expected_in_r;
    end

    % Removal probability
    p_remove = zeros(natom, 1);
    over     = local_density > 1.0;
    p_remove(over) = 1.0 - target_frac ./ local_density(over);
    p_remove = min(p_remove, 1.0 - min_keep);
    p_remove = max(p_remove, 0.0);

    % Remove highest-density atoms first
    [~, order] = sort(local_density, 'descend');
    alive      = true(natom, 1);
    n_removed  = 0;

    for k = 1:natom
        i = order(k);
        if ~alive(i) || p_remove(i) <= 0, continue; end
        if rand() < p_remove(i)
            alive(i) = false;
            n_removed = n_removed + 1;
        end
    end

    % Apply removal
    surviving_ids = atom_ids(alive);
    Atoms = Atoms(alive, :);

    bond_i    = Bonds(:, 2);
    bond_j    = Bonds(:, 3);
    keep_bond = ismember(bond_i, surviving_ids) & ismember(bond_j, surviving_ids);
    n_bonds_rm = sum(~keep_bond);
    Bonds = Bonds(keep_bond, :);

    if ~isempty(Nvec)
        Nvec = Nvec(keep_bond, :);
    end

    obj.log.print(['   [AddDefects:Thin] Removed %d / %d atoms  ' ...
                   '(%d bonds).  %d atoms remain.\n'], ...
        n_removed, natom_before, n_bonds_rm, size(Atoms, 1));
end


%% =======================================================================
%  LOCAL HELPER:  local_bridge_constrictions
%
%  Pass 3 — add single bonds across narrow void constrictions.
%
%  A "constriction" is detected when two atoms sit on opposite sides of
%  a void gap: they are within bridge_max_dist of each other, not already
%  bonded, separated by a region of very low local atom density, and the
%  void extends significantly in the perpendicular direction (confirming
%  the gap is a narrow waist of a larger void rather than just sparse mesh).
%
%  Only the shortest qualifying pair in each local neighbourhood is
%  bridged, preventing over-bridging the same void.
%
%  Reads from obj.defect:
%    bridge_max_dist     - max atom-atom distance to consider for bridging
%                          0 = auto (3.5 × avg node spacing)
%    bridge_void_thresh  - local density at midpoint below this = void
%                          (default 0.25)
%    bridge_perp_width   - min perpendicular void half-width to qualify
%                          0 = auto (1.0 × avg node spacing)
%    bridge_max_degree   - skip atoms already at this bond count
%                          0 = use obj.peratom.Max_peratom_bond
%% =======================================================================
function [Atoms, Bonds, Nvec] = local_bridge_constrictions(obj, Atoms, Bonds, Nvec)

    d = obj.defect;

    rho_atom    = obj.architecture.rho_atom;
    avg_spacing = 1.0 / sqrt(max(rho_atom, 1e-10));
    b           = obj.domain.b;    % fundamental lengthscale

    % All user-facing length settings are in units of b.
    % 0 = auto, computed from node spacing.
    bridge_max_dist = 3.5 * avg_spacing;         % auto default
    if isfield(d, 'bridge_max_dist') && d.bridge_max_dist > 0
        bridge_max_dist = d.bridge_max_dist * b; % user value × b
    end

    void_thresh = 0.25;
    if isfield(d, 'bridge_void_thresh')
        void_thresh = d.bridge_void_thresh;
    end

    perp_width = 1.0 * avg_spacing;              % auto default
    if isfield(d, 'bridge_perp_width') && d.bridge_perp_width > 0
        perp_width = d.bridge_perp_width * b;   % user value × b
    end

    max_deg = obj.peratom.Max_peratom_bond;
    if isfield(d, 'bridge_max_degree') && d.bridge_max_degree > 0
        max_deg = d.bridge_max_degree;
    end

    max_bonds = Inf;   % no cap by default
    if isfield(d, 'bridge_max_bonds') && d.bridge_max_bonds > 0
        max_bonds = d.bridge_max_bonds;
    end

    % Minimum spacing between bridges — controls how densely long voids
    % are bridged. Smaller = more bridges along a snake void.
    % Default: same as bridge_max_dist (conservative — one per neighbourhood).
    bridge_min_spacing = bridge_max_dist;
    if isfield(d, 'bridge_min_spacing') && d.bridge_min_spacing > 0
        bridge_min_spacing = d.bridge_min_spacing * b;
    end

    obj.log.print('   [AddDefects:Bridge] Pass 3 — constriction bridging\n');
    obj.log.print('   [AddDefects:Bridge] max_dist=%.1fb (%.1f)  void_thresh=%.2f  perp_width=%.1fb (%.1f)\n', ...
        bridge_max_dist/b, bridge_max_dist, void_thresh, perp_width/b, perp_width);
    obj.log.print('   [AddDefects:Bridge] max_bonds=%d  min_spacing=%.1fb (%.1f)\n', ...
        max_bonds, bridge_min_spacing/b, bridge_min_spacing);

    atom_ids = Atoms(:, 1);
    % Atoms column layout: [ID | molID | X | Y | Z | deg | nbrs]
    ax       = Atoms(:, 3);
    ay       = Atoms(:, 4);
    natom    = numel(atom_ids);

    % Expected atoms in density check circles
    r_density = 0.8 * avg_spacing;
    expected  = max(rho_atom * pi * r_density^2, 0.5);

    % Compute current degree per atom
    bond_i  = Bonds(:, 2);
    bond_j  = Bonds(:, 3);
    degree  = zeros(natom, 1);
    max_id  = max(atom_ids);
    id2row  = zeros(max_id, 1, 'int32');
    for i = 1:natom
        id2row(atom_ids(i)) = i;
    end
    for b = 1:size(Bonds, 1)
        ri = id2row(bond_i(b));
        rj = id2row(bond_j(b));
        if ri > 0, degree(ri) = degree(ri) + 1; end
        if rj > 0, degree(rj) = degree(rj) + 1; end
    end

    % Build existing bond set for fast lookup (as sorted pair hash)
    existing_pairs = sort([bond_i, bond_j], 2);
    pair_set = containers.Map('KeyType','char','ValueType','logical');
    for b = 1:size(existing_pairs, 1)
        key = sprintf('%d_%d', existing_pairs(b,1), existing_pairs(b,2));
        pair_set(key) = true;
    end

    % Candidate pairs: within bridge_max_dist, not bonded, both have degree room
    new_bonds = zeros(0, size(Bonds, 2));
    new_nvec  = zeros(0, size(Nvec, 2));
    next_bond_id = max(Bonds(:, 1)) + 1;
    default_type = mode(Bonds(:, 5));
    median_N     = median(Nvec(:));

    % For each pair (i,j) we track whether we've already added a bridge
    % in the local neighbourhood — only one per spatial cluster of candidates
    bridged_near = false(natom, 1);   % atoms near a placed bridge

    n_added = 0;

    % Sort candidate atoms by x for a simple sweep
    [~, sorted_i] = sort(ax);

    for si = 1:natom
        i = sorted_i(si);
        if degree(i) >= max_deg, continue; end
        if bridged_near(i),      continue; end

        xi = ax(i);  yi = ay(i);

        % Find nearby atoms within bridge_max_dist
        dx2 = (ax - xi).^2 + (ay - yi).^2;
        candidates = find(dx2 > 0 & dx2 <= bridge_max_dist^2);
        if isempty(candidates), continue; end

        % Sort by distance, try closest first
        [dist_sq, dsort] = sort(dx2(candidates));
        candidates = candidates(dsort);

        for ci = 1:numel(candidates)
            j = candidates(ci);
            if degree(j) >= max_deg, continue; end
            if bridged_near(j),      continue; end

            % Check not already bonded
            id_i = atom_ids(i);
            id_j = atom_ids(j);
            p = sort([id_i, id_j]);
            key = sprintf('%d_%d', p(1), p(2));
            if pair_set.isKey(key), continue; end

            % Check midpoint local density — must be void-like
            mx = (xi + ax(j)) / 2;
            my = (yi + ay(j)) / 2;
            dmid2 = (ax - mx).^2 + (ay - my).^2;
            n_mid = sum(dmid2 <= r_density^2);
            mid_density = n_mid / expected;

            if mid_density >= void_thresh, continue; end

            % Check perpendicular void extent
            % Sample two points offset perpendicular to the bond direction
            dist_ij = sqrt(dist_sq(ci));
            ux = (ax(j) - xi) / dist_ij;   % unit vector i->j
            uy = (ay(j) - yi) / dist_ij;
            % Perpendicular unit vector
            px_u = -uy;  py_u = ux;

            % Check density at perp_width offset in both directions
            p1x = mx + perp_width * px_u;
            p1y = my + perp_width * py_u;
            p2x = mx - perp_width * px_u;
            p2y = my - perp_width * py_u;

            dp1 = (ax - p1x).^2 + (ay - p1y).^2;
            dp2 = (ax - p2x).^2 + (ay - p2y).^2;
            dens_p1 = sum(dp1 <= r_density^2) / expected;
            dens_p2 = sum(dp2 <= r_density^2) / expected;

            % At least one perpendicular direction must also be void-like
            if dens_p1 >= void_thresh && dens_p2 >= void_thresh
                continue;   % void doesn't extend perpendicular — not a constriction
            end

            % Add the bridge bond
            L0 = dist_ij;
            new_row = [next_bond_id, id_i, id_j, L0, default_type];
            new_bonds(end+1, :) = new_row; %#ok<AGROW>
            if ~isempty(Nvec)
                new_nvec(end+1, :) = median_N * ones(1, size(Nvec, 2)); %#ok<AGROW>
            end

            pair_set(key) = true;
            next_bond_id  = next_bond_id + 1;
            degree(i)     = degree(i) + 1;
            degree(j)     = degree(j) + 1;
            n_added       = n_added + 1;

            % Mark neighbourhood of this bridge as already served.
            % Use bridge_min_spacing (not bridge_max_dist) so that
            % multiple bridges can be placed along a long snake void.
            dxj2 = (ax - ax(j)).^2 + (ay - ay(j)).^2;
            near_i = dx2  <= bridge_min_spacing^2;
            near_j = dxj2 <= bridge_min_spacing^2;
            bridged_near = bridged_near | near_i | near_j;

            if n_added >= max_bonds
                obj.log.print('   [AddDefects:Bridge] Cap of %d bonds reached.\n', max_bonds);
                break;
            end

            break;   % one bridge per atom i — move to next
        end
        if n_added >= max_bonds, break; end
    end

    if n_added > 0
        Bonds = [Bonds; new_bonds];
        if ~isempty(Nvec)
            Nvec = [Nvec; new_nvec];
        end
    end

    obj.log.print('   [AddDefects:Bridge] Added %d bridge bonds across constrictions.\n', n_added);
end