function LDpot = ConstructLDPotential(obj, Atoms, Bonds, Nvec)
% -------------------------------------------------------------------------
% ConstructLDPotential
% - Construct local-density potential parameters and table.
%
% Crosslink detection is now molID-based (geometry-agnostic): if the network
% contains more than one unique molID (e.g. bottle-brush with one molID per
% rod), bonds whose two endpoints share a molID are treated as intra-
% molecule (rod) bonds and excluded from the Kuhn-segment sum. Otherwise
% (single-molID networks), every bond contributes. This removes the old
% hardcoded `Bonds(:,5) == 1` check.
%
% Atoms column layout: [ID | molID | X | Y | Z | deg | nbrs...]
%
% INPUT:
%   obj   : network object
%   Atoms : atom array
%   Bonds : bond array
%   Nvec  : per-bond Kuhn segment counts
%
% OUTPUT:
%   LDpot : struct with local-density potential parameters
% -------------------------------------------------------------------------

    %% ------------------------------------------------------------------
    %  0.  Early exit if defects are disabled
    %% ------------------------------------------------------------------
    if ~obj.flags.ipotential
        obj.log.print('   [ConstructLDPotential] Skipped (flags.ipotential = false)\n');
        LDpot = [];
        return;
    end

    if nargin < 4
        error('ConstructLDPotential: requires obj, Atoms, Bonds, and Nvec.');
    end

    Atom_count = size(Atoms,1);
    Bond_count = size(Bonds,1); %#ok<NASGU>

    if Atom_count == 0
        error('ConstructLDPotential: Atoms is empty.');
    end

    if isempty(Nvec)
        error('ConstructLDPotential: Nvec is empty.');
    end

    % ---------------------------------------------------------------------
    % Unpack domain / potential settings
    % ---------------------------------------------------------------------
    xlo = obj.domain.xlo;
    xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;
    yhi = obj.domain.yhi;

    b = obj.domain.b;

    type    = obj.pot.type;
    kLD     = obj.pot.k_LD;
    
    N_rho   = obj.pot.N_rho;
    rho_min = obj.pot.rho_min;
    rho_max = obj.pot.rho_max;

    % Harmonic parameters (none)

    % vdW parameters
    ea      = obj.pot.ea;
    

    if isempty(kLD)
        kLD = 0.414;
    end
    if isempty(N_rho)
        N_rho = 100000;
    end
    if isempty(rho_min)
        rho_min = 0.0;
    end
    if isempty(rho_max)
        rho_max = 500.0;
    end

    if N_rho < 2
        error('ConstructLDPotential: obj.pot.N_rho must be >= 2.');
    end

    drho = (rho_max - rho_min) / (N_rho - 1);
    rho0 = [];

    % ---------------------------------------------------------------------
    % Derived network quantities
    %
    % Crosslink = inter-molecule bond. If the network has multiple molIDs
    % we treat same-molID bonds as intra-molecule (rod-like) and exclude
    % them from the Kuhn-segment sum. For single-molID networks every
    % bond counts.
    % ---------------------------------------------------------------------
    mol_ids = Atoms(:, 2);
    unique_mols = unique(mol_ids(mol_ids > 0));
    multi_mol = numel(unique_mols) > 1;

    if multi_mol && ~isempty(Bonds)
        mol_i = mol_ids(Bonds(:, 2));
        mol_j = mol_ids(Bonds(:, 3));
        crosslink_mask = (mol_i ~= mol_j);
        Total_kuhn_segment = sum(Nvec(crosslink_mask));

        obj.log.print('   LDPot: %d / %d bonds identified as crosslinks (inter-molecule)\n', ...
            sum(crosslink_mask), size(Bonds,1));
    else
        Total_kuhn_segment = sum(Nvec);
    end

    if Total_kuhn_segment <= 0
        error(['ConstructLDPotential: total Kuhn-segment count is <= 0. ' ...
               'This usually means no inter-molecule crosslink bonds remain ' ...
               'after cleanup, or Nvec is zero.']);
    end

    % Determine rod atom size for bottle-brush geometry
    sig_r = 0;
    if strcmpi(obj.architecture.geometry, 'bottle_brush')
        sig_r = obj.architecture.bottlebrush.sigma_c_rod;
    end

    % Calculate sig_c
    sig_c = 0.5 * sqrt(((Total_kuhn_segment / Atom_count) * b^2) + sig_r^2);

    atom_area = (xhi - xlo) * (yhi - ylo);
    if atom_area <= 0
        error('ConstructLDPotential: non-positive in-plane domain area.');
    end

    atom_density = Atom_count / atom_area; %#ok<NASGU>

    % ---------------------------------------------------------------------
    % Construct local-density potential parameters
    % ---------------------------------------------------------------------
 

    if strcmpi(type, 'harmonic')
        R2 = 4.0 * sig_c;
    
        R1 = 0.8 * sig_c;
        rc = 2.0 * sig_c

        rho0 = 0.8 * (R2 / sig_c)^2;
        rho_vec = linspace(rho_min, rho_max, N_rho + 1).';
        pot_density = kLD * (rho_vec - rho0).^2;
    elseif strcmpi(type, 'vdW')
        R2 = 2.85 * sig_c;
    
        R1 = 0.2 * R2;
        rc = 2.0 * sig_c;

        if ea < 0
            error('ConstructLDPotential: vdW potential requires non-negative ea.');
        end

        nu = pi*b^2/4;                                  % Kuhn segment area
        kappa = Total_kuhn_segment / Atom_count;        % Kuhn segments per atom
        Vc = pi*R2^2;                                   % confining volume per atom

        N_bound = Vc / (kappa * nu);
        rho_max = (N_bound - 1) * (1 - 1e-6);
        rho_max
        if rho_min >= N_bound
            error('ConstructLDPotential: vdW rho_min must be below the physical bound %.6g.', N_bound);
        end

        N_points = linspace(rho_min, rho_max, N_rho + 1);
        drho = (rho_max - rho_min) / (N_rho - 1);

        x      = kappa .* (N_points + 1);          % +1 = the central node's own contribution
        term0  = Vc ./ x - nu;                      % Vc/(kappa*(Ntilde+1)) - nu
        term1  = -kLD .* x;
        term2  = ea .* nu .* x.^2 ./ Vc;
        pot_density = term1 .* log(term0) - term2;
    else
        error('ConstructLDPotential: unsupported potential type "%s".', type);
    end

    % ---------------------------------------------------------------------
    % Pack output struct
    % ---------------------------------------------------------------------
    LDpot = struct();

    LDpot.type        = type;
    LDpot.N_LD        = 1;
    LDpot.N_rho       = N_rho;
    LDpot.R_lower     = R1;
    LDpot.R_upper     = R2;
    LDpot.rc          = rc;
    LDpot.rho_min     = rho_min;
    LDpot.rho_max     = rho_max;
    LDpot.drho        = drho;
    LDpot.pot_density = pot_density;
    LDpot.sig_c       = sig_c;

    if strcmpi(type, 'harmonic')
        LDpot.rho0 = rho0;
    end
    if strcmpi(type, 'vdW')
        LDpot.ea = ea;
    end

    obj.log.print('   Constructed LD potential of type %s with %d table entries:\n', ...
        type, numel(pot_density));
    if strcmpi(type, 'harmonic')
        obj.log.print('   Target equilibrium density rho0 = %.4f\n', rho0);
    else
        obj.log.print('   vdW binding energy ea = %.4f\n', ea);
    end
    obj.log.print('   Lower cutoff R1 = %.4f * b\n', R1 / b);
    obj.log.print('   Upper cutoff R2 = %.4f * b\n', R2 / b);
    obj.log.print('   BPM/spring repulsion cutoff rc = %.4f * b\n', rc / b);
end
