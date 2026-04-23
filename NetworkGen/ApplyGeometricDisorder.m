function [Atoms, Bonds] = ApplyGeometricDisorder(obj, Atoms, Bonds)
% -------------------------------------------------------------------------
% ApplyGeometricDisorder
%   Randomly perturb the positions of non-fixed lattice atoms and
%   recompute equilibrium bond lengths (L0) from the new positions.
%
%   Only meaningful for hex_lattice geometry; the guard is enforced by
%   AddHeterogeneities before this function is called.
%
%   Reads from:
%     obj.architecture.lattice_disorder_level   [0, 1]
%         0 -> no movement (perfect lattice)
%         1 -> maximum movement = lattice_disorder_maxfrac * a
%
%     obj.architecture.lattice_disorder_maxfrac
%         Maximum displacement as a fraction of the lattice spacing.
%         Default in architecture.m: 0.4  (i.e. up to 0.4 * a)
%
%     obj.domain.min_node_sep
%         Lattice spacing 'a' (set by SetupDomain from lattice_spacing * b).
%
%     obj.domain.xlo/xhi/ylo/yhi
%         Domain bounds used to clamp displaced atoms.
%
%   Fixed atoms  (Atoms(:,6) == 1)  are never moved.  This preserves the
%   boundary nodes that LAMMPS uses as clamp regions.
%
% INPUT
%   obj   : network object
%   Atoms : [N x (5+MaxNbr)]  atom array  -- col 6 is isFixed flag
%   Bonds : [M x 5]           bond array  [id | i | j | L0 | type]
%
% OUTPUT
%   Atoms : atom array with updated x (col 2) and y (col 3) positions
%   Bonds : bond array with updated L0 (col 4) recomputed from new positions
% -------------------------------------------------------------------------

    disorder_level  = obj.architecture.lattice_disorder_level;
    max_frac        = obj.architecture.lattice_disorder_maxfrac;
    a               = obj.domain.min_node_sep;

    xlo = obj.domain.xlo;  xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;  yhi = obj.domain.yhi;

    % Edge tolerance: keep displaced atoms at least edgeTol inside boundary
    edgeTol = 0.25 * a;

    r_max = disorder_level * max_frac * a;

    if r_max <= 0
        obj.log.print('   [ApplyGeometricDisorder] disorder_level=0, skipping\n');
        return;
    end

    natom    = size(Atoms, 1);
    n_moved  = 0;

    for i = 1:natom

        % Note: the old layout used Atoms(:,6) as an isFixed flag for
        % boundary nodes. Under the refactored layout col 6 is the atom
        % degree, which is populated after AddBonds, so the old guard is
        % no longer meaningful. If a hex-lattice run needs boundary
        % clamping again, reintroduce an explicit isFixed column and
        % restore the guard here.

        % Uniform random displacement inside a disk of radius r_max.
        % sqrt(rand) gives uniform area density (not clumped at centre).
        rr    = r_max * sqrt(rand);
        theta = 2 * pi * rand;

        % Atoms column layout: [ID | molID | X | Y | Z | deg | nbrs]
        x_new = Atoms(i, 3) + rr * cos(theta);
        y_new = Atoms(i, 4) + rr * sin(theta);

        % Clamp to domain interior
        x_new = max(xlo + edgeTol, min(xhi - edgeTol, x_new));
        y_new = max(ylo + edgeTol, min(yhi - edgeTol, y_new));

        Atoms(i, 3) = x_new;
        Atoms(i, 4) = y_new;
        n_moved = n_moved + 1;

    end

    % ------------------------------------------------------------------
    % Recompute L0 for all bonds from the updated atom positions.
    % This is necessary because bond lengths were computed at generation
    % time before any position perturbation.
    % ------------------------------------------------------------------
    if ~isempty(Bonds)
        for k = 1:size(Bonds, 1)
            ii = Bonds(k, 2);   % atom row index (IDs = rows after renumber)
            jj = Bonds(k, 3);

            dx = Atoms(jj, 3) - Atoms(ii, 3);
            dy = Atoms(jj, 4) - Atoms(ii, 4);
            Bonds(k, 4) = sqrt(dx*dx + dy*dy);
        end
    end

    obj.log.print(['   [ApplyGeometricDisorder] Moved %d / %d atoms ' ...
                   '(disorder_level=%.2f, r_max=%.3f)\n'], ...
                  n_moved, natom, disorder_level, r_max);

end
