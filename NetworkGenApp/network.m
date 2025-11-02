classdef network < handle

properties
    % general options
    dist_type
    Nreplicates
    b
    Lx
    Ly
    imanualseed
    seed
    iplot
    isave
    lammps_data_file
    bond_table_file
    write_location
    iadvancedoptions
    % polydisperse options
    distribution_assignment_mode_poly

    % bimodal options
    N1
    N2
    bin_window_method
    distribution_assignment_mode
    distribution_height_mode
    long_first
    P
    N2_bonds
    dr1
    dr2
    
end

properties (SetObservable)
    % observable log property: changes fire PostSet which we forward to
    % the logUpdated event via a small constructor listener below
    log
end

events
    logUpdated
end

methods

    function obj = network()
        % Constructor: forward changes to the observable `log` property as
        % the `logUpdated` event so external listeners can subscribe to
        % either the property's PostSet or the event.
        addlistener(obj,'log','PostSet',@(src,ev) notify(obj,'logUpdated'));
    end

    function [Domain,Atoms,Bonds,Nvec] = generateNetwork(obj)

        % prepare options structure
        options.dist_type          = obj.dist_type;
        options.Nreplicates        = obj.Nreplicates;
        options.b                  = obj.b;
        options.Lx                 = obj.Lx;
        options.Ly                 = obj.Ly;
        options.imanualseed        = obj.imanualseed;
        options.seed               = obj.seed;
        options.iplot              = obj.iplot;
        options.isave              = obj.isave;
        options.lammps_data_file   = obj.lammps_data_file;
        options.bond_table_file    = obj.bond_table_file;
        options.write_location     = obj.write_location;

        % A. Polydisperse options
        % ------------------------------------------------------------------
        % --- mode selection ---
        options.polydisperse.distribution_assignment_mode = obj.distribution_assignment_mode_poly;   % 'geom' | 'range' | 'pmf'

        % --- shared / guards ---
        options.polydisperse.min_N           = 1;        % lower bound for all modes
        options.polydisperse.align_to_length = 'ascend'; % 'ascend' (shortest→smallest N) | 'none'
        options.polydisperse.kuhn_rounding   = 'round';  % 'round' | 'ceil' | 'floor' (used in 'geom')

        % --- 'geom' mode (N ≈ (L/b)^2) ---
        % uses: b (global), kuhn_rounding, min_N

        % --- 'range' mode (map lengths → [N_min,N_max]) ---
        options.polydisperse.N_range_method  = 'rank';   % 'rank' | 'linear'
        options.polydisperse.N_target_min    = 20;       % integer lower target
        options.polydisperse.N_target_max    = 120;      % integer upper target

        % --- 'pmf' mode (truncated geometric with hard cap based on exp distribution) ---
        options.polydisperse.pmf_nu0         = 20;       % ν0 (minimum)
        options.polydisperse.pmf_meanN       = 30;       % target mean of ν after truncation
        options.polydisperse.pmf_cut_mode    = 'cap';    % keep as 'cap'
        options.polydisperse.pmf_nu_max      = 60;       % hard maximum ν (≥ ν0)
        options.polydisperse.integerize_rule = 'largest_remainder'; % allocation method

        % B. Bimodal options
        % ------------------------------------------------------------------
        options.bimodal.N1                 = obj.N1;
        options.bimodal.N2                 = obj.N2;
        options.bimodal.dr1                = obj.dr1;
        options.bimodal.dr2                = obj.dr2;
        options.bimodal.min_N              = 1;        % lower bound for all modes
        options.bimodal.kuhn_rounding      = 'round';  % 'round' | 'ceil' | 'floor' (used in 'geom')

        options.bimodal.long_first         = obj.long_first;    
        options.bimodal.bin_window_method  = obj.bin_window_method; 

        % --- mode selection ---
        % 'single' mode: applies N1 and N2 directly based on geometry (N ≈ (L/b)^2)
        % 'range' mode: uses: b (global), kuhn_rounding, min_N
        options.bimodal.distribution_assignment_mode = obj.distribution_assignment_mode; % 'single' | 'geom' | 'gaussian'

        % --- distribution height mode selection ---
        options.bimodal.distribution_height_mode = obj.distribution_height_mode; % 'prob' | 'fixed'

        % --- 'prob' mode ---
        options.bimodal.P = obj.P;          % desired fraction of type 2 bonds

        % --- 'fixed' mode ---
        options.bimodal.N2_number = obj.N2_bonds;  % fraction of type 2 bonds

        % C. Additional advanced options
        % ------------------------------------------------------------------
        advancedOptions.iadvancedoptions = obj.iadvancedoptions;
        if advancedOptions.iadvancedoptions
            advancedOptions.Rho_atom = 0.0078;
            advancedOptions.Max_peratom_bond = 5;
            advancedOptions.bond_global_try_limit_multiplier = 200;
            advancedOptions.max_attempts_without_progress_multiplier = 10;
            advancedOptions.min_degree_keep = 1;
            advancedOptions.cutoff_multiply = 6; %units of b

            %Bimodal advanced
            %advancedOptions.bimodal.bin_std1_factor = 0.4;
            %advancedOptions.bimodal.bin_std2_factor = 0.15;
            %advancedOptions.bimodal.bin_width1_factor = 2.355;
            %advancedOptions.bimodal.bin_width2_factor = 2.355;
        end

        obj.log = sprintf("Starting network generation...\n");

        for ii = 1:options.Nreplicates;

            newline = sprintf('----------------------------------------\n');
            obj.log = append(obj.log, newline);

            fprintf('Generating network replicate %d of %d...\n',ii,options.Nreplicates);
            newline = sprintf('Generating network replicate %d of %d...\n',ii,options.Nreplicates);
            obj.log = append(obj.log, newline);

            %% A. Prepare replicate-specific information
            % 1. Set seed
            if options.imanualseed
                if length(options.seed) ~= options.Nreplicates
                    newline = sprintf('Error: not enough manual seeds provided for the number of replicates Expected %d, got %d\n', options.Nreplicates, length(options.seed));
                    obj.log = append(obj.log, newline);
                    error('Error: not enough manual seeds provided for the number of replicates Expected %d, got %d', options.Nreplicates, length(options.seed));
                end
                options.seed = options.seed(ii);
            else
                options.seed = randi(1e6);
            end
            rng(options.seed);

            % 2. Set replicate-specific file names
            if options.Nreplicates > 1
                options.lammps_data_file   = sprintf('%04d_%s',ii,obj.lammps_data_file);
                %options.lammps_visual_file = sprintf('%04d_%s',ii,options.lammps_visual_file);
                options.bond_table_file    = sprintf('%04d_%s',ii,obj.bond_table_file);
            end

            %% B. Setup the system
            [Domain,obj] = NetworkGenSetup(options,advancedOptions,obj);

            %% C. Scatter nodes with minimum spacing
            [Atoms,obj] = NetworkGenScatterNodes(Domain,obj);

            %% D. Connect nodes within bonds
            if strcmp(options.dist_type,'polydisperse')
                % Polydisperse network
                [Atoms,Bonds,obj] = NetworkGenConnectNodesPolydisperse(Domain,Atoms,options,obj);
                [Nvec,obj] = NetworkGenAssignKuhnPolydisperse(Bonds, options,obj);
            elseif strcmp(options.dist_type,'bimodal')
                % Bimodal network
                [Atoms,Bonds,obj] = NetworkGenConnectNodesBimodal(Domain,Atoms,options,advancedOptions,obj);
                [Nvec] = NetworkGenAssignKuhnBimodal(Bonds, options,obj);
            else
                newline = sprintf('Error: distribution type: %s not recognized\n', options.dist_type);
                obj.log = append(obj.log, newline);
                error('Error: distribution type: %s not recognized', options.dist_type);
            end

            %% E. Save data files
            if options.isave
                [obj] = NetworkGenWriteDataFiles(Domain,Atoms,Bonds,Nvec,options,obj);
            else
                newline = sprintf('   Data files not written because options.isave = false\n');
                obj.log = append(obj.log, newline);
            end

        end

        newline = sprintf('----------------------------------------\n\n');
        obj.log = append(obj.log, newline);

        newline = sprintf("Done!\n\n");
        obj.log = append(obj.log, newline);

        newline = sprintf('Visualizing last generated network... (this may take a minute)\n');
        obj.log = append(obj.log, newline);
    end

end
end