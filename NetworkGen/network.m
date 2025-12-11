classdef network < handle

properties

    % Define properites and set default parameters
    %%% --------------------------- Global settings ----------------------
    network_geometry    = 'random';     % 'random' | 'hex_lattice'
    dist_type           = 'bimodal';    % 'bimodal' | 'polydisperse'
    Nreplicates         = 1;            % number of network replicates to generate
    b                   = 1.6;          % Kuhn length
    Lx                  = 150*2;        % Domain size in x (in units of b)
    Ly                  = 90*2;         % Domain size in y (in units of b)
    scale               = 1.8;          % scale factor for final domain size
    boundary_box        = 'fixed';      % 'fixed' | 'periodic'
    imanualseed         = false;        % 'false' | 'true' use manual seeds
    seed                = [42];         % manual seed values
    iplot               = true;         % 'true' | 'false' show visualization
    isave               = true;         % 'true' | 'false' save data files
    write_location      = './networks'; % location to write output files
    lammps_data_file    = 'PolyNetwork';% prefix file name for LAMMPS data output
    lammps_visual_file  = 'PolyVisual'; % prefix file name for LAMMPS visualization output
    bond_table_file     = 'bond';       % file name for bond table output
    save_name_mode      = true;         % 'true' | 'false' auto naming mode for output files
    smp_number          = 1;            % sample number for naming output files

    % LD parameters
    kLD                 = 2.0*4.14; % strength factor
    N_rho               = 100000;   % number of density points
    rho_min             = 0.0;      % minimum density
    rho_max             = 500;      % maximum density

    % Lattice options
    lattice_a                   = 6;    % lattice spacing (in units of b)
    lattice_disorder_level      = 0.1;  % fraction of lattice spacing for positional disorder
    lattice_disorder_max_frac_a = 0.4;  % max displacement radius as fraction of 'a'
    lattice_topo_disorder_flag  = true; % 'true' | 'false' enable/disable bond deletions
    lattice_max_del_per_node    = 1;    % max number of bonds to delete per node
    lattice_min_degree_keep     = 5;    % don't let any node go below this degree

    % polydisperse options
    distribution_assignment_mode_poly   = 'mono';       % 'geom' | 'range' | 'pmf' | 'mono'
    min_N_poly                          = 1;            % lower bound for all modes
    align_to_length_poly                = 'ascend';     % 'ascend' (shortest→smallest N) | 'none'
    kuhn_rounding_poly                  = 'round';      % 'round' | 'ceil' | 'floor' (used in 'geom')
    N_range_method_poly                 = 'rank';       % 'rank' | 'linear'
    N_target_min_poly                   = 20;           % integer lower target
    N_target_max_poly                   = 120;          % integer upper target
    pmf_nu0_poly                        = 50;           % ν0 (minimum)
    pmf_meanN_poly                      = 70;           % target mean of ν after truncation
    pmf_cut_mode_poly                   = 'cap';        % keep as 'cap'
    pmf_nu_max_poly                     = 280;          % hard maximum ν (≥ ν0)
    integerize_rule_poly                = 'largest_remainder'; % allocation method

    % bimodal options
    N1 = 35; % Number of Kuhn segments in type 1 bonds
    N2 = 60; % Number of Kuhn segments in type 2 bonds
    
    bin_window_method                    = 'manual';    % Method for determining bin width of bimodal dist: 'manual' or 'adaptive'
    manual_deviation_type                = 'mixed';     % For manual bins is the standard deviation of bin width: 'kuhn' or 'mixed'
    distribution_assignment_mode_bimodal = 'gaussian';  % Kuhn segment assigment method: 'single' or 'geom' or 'gaussian'
    distribution_height_mode             = 'prob';      % Distribution height method: 'prob' or 'fixed'
    
    long_first          = true;         % 'true' | 'false' enable long-first mode
    double_network_flag = true;         % 'true' | 'false' enable double network style
    auto_N2_flag        = true;         % 'true' | 'false' automatically overrides N2 given the spacing ratio and desired pre-stretch
    alpha               = 1.0;          % spacing ratio of large mesh to small mesh
    P                   = 1.0;          % Prob: desired fraction of type 2 bonds
    N2_bonds            = 2;            % Fixed: desired number of type 2 bonds
    lam1                = 0.17;         % Prestretched length of type 1 bonds: lam1 = [0 1], 1/sqrt(N1) (default)
    lam2                = 0.17;         % Prestretched length of type 2 bonds: lam2 = [0 1], 1/sqrt(N2) (default)
    stdN1               = 10;           % std of N1 Kuhn segment distribution  
    stdN2               = 10;           % std of N2 Kuhn segment distribution 
    stdr1               = 3;            % std of the end-to-end length for r1;
    stdr2               = 10;           % std of the end-to-end length for r2;

    min_N_bimodal           = 1;        % lower bound for all modes
    kuhn_rounding_bimodal   = 'round';  % 'round' | 'ceil' | 'floor' (used in 'geom')

    % advanced options
    iadvancedoptions    = false;
    rho_atom            = 0.0078;
    max_peratom_bond    = 5;
    bond_global_try_limit_multiplier            = 200;
    max_attempts_without_progress_multiplier    = 10;
    min_degree_keep     = 1;
    cutoff_multiply     = 6; %units of b
end

methods

    function [Domain,Atoms,Bonds,Nvec,order] = generateNetwork(obj)

        %%% Options structure is now redundant with object properties; still keep for compatibility with old code
        % prepare options structure

        % prepare options structure
        options.dist_type          = obj.dist_type;
        options.network_geometry   = obj.network_geometry;
        options.Nreplicates        = obj.Nreplicates;
        options.b                  = obj.b;
        options.Lx                 = obj.Lx;
        options.Ly                 = obj.Ly;
        options.boundary_box       = obj.boundary_box;
        options.imanualseed        = obj.imanualseed;
        options.seed               = obj.seed;
        options.iplot              = obj.iplot;
        options.isave              = obj.isave;
        options.lammps_data_file   = obj.lammps_data_file;
        options.lammps_visual_file = obj.lammps_visual_file;
        options.bond_table_file    = obj.bond_table_file;
        options.write_location     = obj.write_location;
        options.save_name_mode     = obj.save_name_mode;
        options.smp_number         = obj.smp_number;

        % Local density pot options
        options.LDpot_strength     = obj.kLD;        % strength factor
        options.LDpot_N_rho        = obj.N_rho;      % number of density points
        options.LDpot_rho_min      = obj.rho_min;    % minimum density
        options.LDpot_rho_max      = obj.rho_max;    % maximum density

        % Lattice-specific options
        options.lattice.a                  = obj.lattice_a*obj.b;
        options.lattice.edgeTol            = 0.25*obj.lattice_a;
        options.lattice.disorder_level     = obj.lattice_disorder_level;
        options.lattice.disorder_max_frac_a= obj.lattice_disorder_max_frac_a;

        options.lattice.enable_topo_disorder   = obj.lattice_topo_disorder_flag;
        options.lattice.max_topo_del_per_node  = obj.lattice_max_del_per_node;
        options.lattice.min_degree_keep        = obj.lattice_min_degree_keep;

        % Distribution options
        % A. Polydisperse options
        % ------------------------------------------------------------------
        % --- mode selection ---
        options.polydisperse.distribution_assignment_mode = obj.distribution_assignment_mode_poly;   % 'geom' | 'range' | 'pmf'

        % --- shared / guards ---
        options.polydisperse.min_N           = obj.min_N_poly;        
        options.polydisperse.align_to_length = obj.align_to_length_poly;  % 'ascend' (shortest→smallest N) | 'none'
        options.polydisperse.kuhn_rounding   = obj.kuhn_rounding_poly;    % 'round' | 'ceil' | 'floor' (used in 'geom')

        % --- 'geom' mode (N ≈ (L/b)^2) ---
        % uses: b (global), kuhn_rounding, min_N

        % --- 'range' mode (map lengths → [N_min,N_max]) ---
        options.polydisperse.N_range_method  = obj.N_range_method_poly;   % 'rank' | 'linear'
        options.polydisperse.N_target_min    = obj.N_target_min_poly;     % integer lower target
        options.polydisperse.N_target_max    = obj.N_target_max_poly;     % integer upper target

        % --- 'pmf' mode (truncated geometric with hard cap based on exp distribution) ---
        options.polydisperse.pmf_nu0         = obj.pmf_nu0_poly;          % ν0 (minimum)
        options.polydisperse.pmf_meanN       = obj.pmf_meanN_poly;        % target mean of ν after truncation
        options.polydisperse.pmf_cut_mode    = obj.pmf_cut_mode_poly;     % keep as 'cap'
        options.polydisperse.pmf_nu_max      = obj.pmf_nu_max_poly;       % hard maximum ν (≥ ν0)
        options.polydisperse.integerize_rule = obj.integerize_rule_poly;   % allocation method

        % B. Bimodal options
        % ------------------------------------------------------------------
        options.bimodal.N1                 = obj.N1;
        options.bimodal.N2                 = obj.N2;
        options.bimodal.std1               = obj.stdN1;
        options.bimodal.std2               = obj.stdN2;
        options.bimodal.stdr1              = obj.stdr1;
        options.bimodal.stdr2              = obj.stdr2;
        options.bimodal.lam1               = obj.lam1;
        options.bimodal.lam2               = obj.lam2;
        options.bimodal.min_N              = obj.min_N_bimodal;          % lower bound for all modes
        options.bimodal.kuhn_rounding      = obj.kuhn_rounding_bimodal;  % 'round' | 'ceil' | 'floor' (used in 'geom')

        options.bimodal.long_first         = obj.long_first;    
        options.bimodal.bin_window_method  = obj.bin_window_method; 
        options.bimodal.deviation_type     = obj.manual_deviation_type;
        options.double_network.flag        = obj.double_network_flag;
        options.double_network.autoN2      = obj.auto_N2_flag;
        options.double_network.alpha       = obj.alpha;

        % --- mode selection ---
        % 'single' mode: applies N1 and N2 directly based on geometry (N ≈ (L/b)^2)
        % 'range' mode: uses: b (global), kuhn_rounding, min_N
        options.bimodal.distribution_assignment_mode = obj.distribution_assignment_mode_bimodal; % 'single' | 'geom' | 'gaussian'

        % --- distribution height mode selection ---
        options.bimodal.distribution_height_mode = obj.distribution_height_mode; % 'prob' | 'fixed'

        % --- 'prob' mode ---
        options.bimodal.P = obj.P;          % desired fraction of type 2 bonds

        % --- 'fixed' mode ---
        options.bimodal.N2_number = obj.N2_bonds;  % fraction of type 2 bonds

        % C. Additional advanced options
        % ------------------------------------------------------------------
        advancedOptions.iadvancedoptions = obj.iadvancedoptions;
        if obj.iadvancedoptions
            advancedOptions.Rho_atom = obj.rho_atom;
            advancedOptions.Max_peratom_bond = obj.max_peratom_bond;
            advancedOptions.bond_global_try_limit_multiplier = obj.bond_global_try_limit_multiplier;
            advancedOptions.max_attempts_without_progress_multiplier =obj.max_attempts_without_progress_multiplier;
            advancedOptions.min_degree_keep = obj.min_degree_keep;
            advancedOptions.cutoff_multiply = obj.cutoff_multiply; %units of b

        end

        %% Loop over replicates
        for ii = 1:obj.Nreplicates
            fprintf('Generating network replicate %d of %d...\n',ii,obj.Nreplicates);
    
            %% A. Prepare replicate-specific information
            % 1. Set seed
            if obj.imanualseed
                if length(obj.seed) ~= obj.Nreplicates
                    error('Error: not enough manual seeds provided for the number of replicates Expected %d, got %d', obj.Nreplicates, length(obj.seed));
                end
                options.seed = obj.seed(ii);
            else
                options.seed = randi(1e6);
            end
            rng(options.seed);

            % 2. Set sample name
            sample_suffix = sprintf('SMP%04d',obj.smp_number);
            options.sample_suffix = sample_suffix;
    
            % 3. Set replicate-specific file names
            if obj.Nreplicates > 1
                replicate_suffix = sprintf('N%04d_',ii);
                options.replicate_suffix = replicate_suffix;
            else
                replicate_suffix = sprintf('N%04d',1);
                options.replicate_suffix = replicate_suffix;
            end

            %% B. Setup the system
            [Domain] = NetworkGenSetup(options,advancedOptions);

            %% C. Add crosslink nodes and connect with bonds
            %1. Check geometry type
            if strcmp(options.network_geometry, 'random')

                % ---- Random Network Generation ----
                % Add atoms
                [Atoms] = NetworkGenScatterNodes(Domain);
        
                % Add bonds
                if strcmp(obj.dist_type,'polydisperse')
                    [Atoms,Bonds] = NetworkGenConnectNodesPolydisperse(Domain,Atoms,options);
                    [Nvec]        = NetworkGenAssignKuhnPolydisperse(Bonds,Atoms, options);
                elseif strcmp(obj.dist_type,'bimodal')
                    [Atoms,Bonds,options] = NetworkGenConnectNodesBimodal(Domain,Atoms,options,advancedOptions);
                    [Nvec]                = NetworkGenAssignKuhnBimodal(Bonds,Atoms, options);
                else
                    error('Error: distribution type: %s not recognized', obj.dist_type);
                end

            elseif strcmp(options.network_geometry, 'hex_lattice')

                % ---- Hexagonal lattice ----
                % Add atoms
                [Atoms, latticeData] = NetworkGenLatticeScatterNodes(Domain, options);

                % Add bonds
                [Atoms, Bonds] = NetworkGenLatticeConnectNodes(Domain, Atoms, latticeData, options);
        
                % Add disorder
                if isfield(options,'lattice') && isfield(options.lattice,'enable_topo_disorder')
                    if options.lattice.enable_topo_disorder
                        [Atoms, Bonds] = NetworkApplyLatticeTopoDisorder(Atoms, Bonds, options);
                    end
                end

                % Assign Kuhn segments
                if strcmp(obj.dist_type,'polydisperse')
                    [Nvec] = NetworkGenAssignKuhnPolydisperse(Bonds,Atoms, options);
                elseif strcmp(obj.dist_type,'bimodal')
                    [Nvec] = NetworkGenAssignKuhnBimodal(Bonds,Atoms, options);
                else
                    error('Error: distribution type: %s not recognized', obj.dist_type);
                end

            else
                error('Unknown network_geometry: %s', options.network_geometry);
            end

            %2. Contruct local density potential
            [LDpot] = NetworkGenConstructLDPotential(Domain, Atoms, Bonds, Nvec, options);
    
            %% D. Scale domain if needed
            if (obj.scale ~= 1.0)
                [Domain, Atoms, Bonds] = NetworkScaleDomain(Domain, Atoms, Bonds, obj.scale);
            end

            %% E. Show visualization and statistics
            NetworkGenVisualize(Domain, Atoms, Bonds, Nvec, obj.scale, options);
    
            %% F. Compute order parameters
            [order] = NetworkComputeOrder(Atoms,Bonds);
    
            %% G. Write data files
            NetworkGenWriteDataFiles(Domain, Atoms, Bonds, Nvec, LDpot, options, order);
    
        end

        fprintf('Done!\n')
        
    end

end
end