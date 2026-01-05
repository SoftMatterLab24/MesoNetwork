% -------------------------------------------------------------------------
% Mesoscale network generator for polydisperse and bimodal networks
%  (2D, no PBC) — MATLAB R2016a
% - Enforces minimum spacing between scattered nodes (rejection sampling)
% - Randomly connects nearby nodes under distance cutoff (with guards)
% - Assigns bonds distribution based on desired statistics
% - Iteratively prunes nodes with degree <= 1
% - Exports:
%     * Network.txt  (LAMMPS data: atoms/bonds; atom type=1, bond type=1)
%     * bond.table   (# Chain stats; KEY; N <#bonds>; lines: id i j N b)
% - Per-bond Kuhn segments specified: N1, N2 
% -------------------------------------------------------------------------

clc; clear; close all;
warning off backtrace  % disable stack trace for warnings

%% --------------------------- Global settings ----------------------
% --- Network geometry: 'random' or 'hex_lattice'
network_geometry = 'hex_lattice';   % 'random' | 'hex_lattice'

% --- Distribution type: 'bimodal' or 'polydisperse'
dist_type = 'polydisperse';

% --- Number of networks to generate
Nreplicates = 1;

% --- Kuhn length
b = 1.6;          % Kuhn length (in nm)

% --- Domain size
Lx = 150*8;       % Domain size in x (in units of b)
Ly = 90*8;        % Domain size in y (in units of b)

% --- Domain size scaler
scale = 3.4;  % e.g., halve the system dimensions

% --- Boundary Conditions
boundary_box = 'fixed'; % 'fixed' or 'periodic' boundaries

% --- Seed options
imanualseed = false;  % true: manual seed; false: random seed
seed = [1];

% --- Visualization
iplot = true;    % Show 

% --- Save options
isave              = true;                                    % Save data files
write_location     = './networks';                            % Location to write output files
lammps_data_file   = 'PolyNetwork';                           % Prefix file name for LAMMPS data output
lammps_visual_file = 'PolyVisual';                            % Prefix file name for LAMMPS visualization output
bond_table_file    = 'bond';                                  % File name for bond table output   
save_name_mode     = true;                                    % true: auto add sample info to file names; false: use only fixed names
smp_number         = 1;                                       % Sample number for file naming <- to be used for auto input script making (data sweeps)

%% --------------------- Local Density Potential Options  -----------------
kLD     = 2*4.14; % strength factor
N_rho   = 100000; % number of density points
rho_min = 0.0;    % minimum density
rho_max = 500;    % maximum density

%% --------------------- Hexagonal Lattice Options ------------------------
% --- Lattice spacing (center-to-center distance between neighboring nodes)
lattice_a = 6 * b;   % you can tune this

% --- Lattice disorder (0 = perfect geometry, 1 = strong geometric disorder)
lattice_disorder_level      = 1;   % try 0, 0.3, 0.6, 1.0
lattice_disorder_max_frac_a = 0.4;  % max displacement radius as fraction of 'a'

% --- Topological disorder options (bond deletions)
lattice_topo_disorder_flag      = true;  % enable/disable bond deletions
lattice_max_del_per_node        = 1;     % max bonds to attempt to delete per node at disorder=1
lattice_min_degree_keep         = 5;     % don't let any node go below this degree

%%% --------------------- Distribution OPTIONS -------------------------%%%
%% --------------------- Polydisperse Options -----------------------------

distribution_assignment_mode_poly = 'mono';  % Kuhn segment assigment method: 'geom' | 'range' | 'pmf' | 'mono'

%% --------------------- Bimodal Options ---------------------------
% --- Average chain Kuhn segments
N1 = 35; 
N2 = 60; 

% --- 
bin_window_method            = 'manual';    % Method for determining bin width of bimodal dist: 'manual' or 'adaptive'
manual_deviation_type        = 'mixed';     % For manual bins is the standard deviation of bin width: 'kuhn' or 'mixed'
distribution_assignment_mode = 'gaussian';  % Kuhn segment assigment method: 'single' or 'geom' or 'gaussian'
distribution_height_mode     = 'prob';      % Distribution height method: 'prob' or 'fixed'
long_first                   =  true;       % enable long-first mode

% --- Double network params
double_network_flag = true;                 % enable double network style
auto_N2_flag        = false;                % automatically overrides N2 given the spacing ratio and desired pre-stretch
alpha = 2.5;                                % spacing ratio of large mesh to small mesh

% --- Height mode settings (only one is used)
P = 1.0;      % Prob: desired fraction of type 2 bonds
N2_bonds = 2; % Fixed: desired number of type 2 bonds

%%% Manual mode settings (only used if bin_window_method = manual)
% --- Prestretch
lam1 = 1/sqrt(N1);   % Prestretched length of type 1 bonds: lam1 = [0 1], 1/sqrt(N1) (default)
lam2 = 0.35;         % Prestretched length of type 2 bonds: lam2 = [0 1], 1/sqrt(N2) (default)

% NOTE: Kuhn uses only kuhn, mixed uses both
% --- Deviation in Kuhn segment
stdN1 = 10; % std of N1 Kuhn segment distribution         
stdN2 = 10; % std of N2 Kuhn segment distribution 

% --- Deviation in end-to-end length
stdr1 = 2;   % std of the end-to-end length for r1;
stdr2 = 12;  % std of the end-to-end length for r2; [2 5 15 25]

%% --------------------- Advanced Options --------------------------
iadvancedoptions = false;

%% --------------------- Network Generation ------------------------
%%% -----------  !!!DO NOT EDIT BELOW THIS LINE!!! ----------- %%%

% prepare options structure
options.dist_type          = dist_type;
options.network_geometry   = network_geometry;
options.Nreplicates        = Nreplicates;
options.b                  = b;
options.Lx                 = Lx;
options.Ly                 = Ly;
options.boundary_box       = boundary_box;
options.imanualseed        = imanualseed;
options.seed               = seed;
options.iplot              = iplot;
options.isave              = isave;
options.lammps_data_file   = lammps_data_file;
options.lammps_visual_file = lammps_visual_file;
options.bond_table_file    = bond_table_file;
options.write_location     = write_location;
options.save_name_mode     = save_name_mode;
options.smp_number         = smp_number;

% Local density pot options
options.LDpot_strength     = kLD;        % strength factor
options.LDpot_N_rho        = N_rho;      % number of density points
options.LDpot_rho_min      = rho_min;    % minimum density
options.LDpot_rho_max      = rho_max;    % maximum density

% Lattice-specific options
options.lattice.a                  = lattice_a;
options.lattice.edgeTol            = 0.25*lattice_a;
options.lattice.disorder_level     = lattice_disorder_level;
options.lattice.disorder_max_frac_a= lattice_disorder_max_frac_a;

options.lattice.enable_topo_disorder   = lattice_topo_disorder_flag;
options.lattice.max_topo_del_per_node  = lattice_max_del_per_node;
options.lattice.min_degree_keep        = lattice_min_degree_keep;

% Distribution options
% A. Polydisperse options
% ------------------------------------------------------------------
% --- mode selection ---
options.polydisperse.distribution_assignment_mode = distribution_assignment_mode_poly;   % 'geom' | 'range' | 'pmf'

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
options.polydisperse.pmf_nu0         = 40;       % ν0 (minimum)
options.polydisperse.pmf_meanN       = 50;       % target mean of ν after truncation
options.polydisperse.pmf_cut_mode    = 'cap';    % keep as 'cap'
options.polydisperse.pmf_nu_max      = 280;       % hard maximum ν (≥ ν0)
options.polydisperse.integerize_rule = 'largest_remainder'; % allocation method

% B. Bimodal options
% ------------------------------------------------------------------
options.bimodal.N1                 = N1;
options.bimodal.N2                 = N2;
options.bimodal.std1               = stdN1;
options.bimodal.std2               = stdN2;
options.bimodal.stdr1              = stdr1;
options.bimodal.stdr2              = stdr2;
options.bimodal.lam1               = lam1;
options.bimodal.lam2               = lam2;
options.bimodal.min_N              = 1;        % lower bound for all modes
options.bimodal.kuhn_rounding      = 'round';  % 'round' | 'ceil' | 'floor' (used in 'geom')

options.bimodal.long_first         = long_first;    
options.bimodal.bin_window_method  = bin_window_method; 
options.bimodal.deviation_type     = manual_deviation_type;
options.double_network.flag        = double_network_flag;
options.double_network.autoN2      = auto_N2_flag;
options.double_network.alpha       = alpha;

% --- mode selection ---
% 'single' mode: applies N1 and N2 directly based on geometry (N ≈ (L/b)^2)
% 'range' mode: uses: b (global), kuhn_rounding, min_N
options.bimodal.distribution_assignment_mode = distribution_assignment_mode; % 'single' | 'geom' | 'gaussian'

% --- distribution height mode selection ---
options.bimodal.distribution_height_mode = distribution_height_mode; % 'prob' | 'fixed'

% --- 'prob' mode ---
options.bimodal.P = P;          % desired fraction of type 2 bonds

% --- 'fixed' mode ---
options.bimodal.N2_number = N2_bonds;  % fraction of type 2 bonds

% C. Additional advanced options
% ------------------------------------------------------------------
advancedOptions.iadvancedoptions = iadvancedoptions;
if iadvancedoptions
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

%% Loop over replicates
for ii = 1:Nreplicates
    fprintf('Generating network replicate %d of %d...\n',ii,Nreplicates);
    
    %% A. Prepare replicate-specific information
    % 1. Set seed
    if imanualseed
        if length(seed) ~= Nreplicates
            error('Error: not enough manual seeds provided for the number of replicates Expected %d, got %d', Nreplicates, length(seed));
        end
        options.seed = seed(ii);
    else
        options.seed = randi(1e6);
    end
    rng(options.seed);

    % 2. Set sample name
    sample_suffix = sprintf('SMP%04d',smp_number);
    options.sample_suffix = sample_suffix;
    
    % 3. Set replicate-specific file names
    if Nreplicates > 1
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
        if strcmp(dist_type,'polydisperse')
            [Atoms,Bonds] = NetworkGenConnectNodesPolydisperse(Domain,Atoms,options);
            [Nvec]        = NetworkGenAssignKuhnPolydisperse(Bonds,Atoms, options);
        elseif strcmp(dist_type,'bimodal')
            [Atoms,Bonds,options] = NetworkGenConnectNodesBimodal(Domain,Atoms,options,advancedOptions);
            [Nvec]                = NetworkGenAssignKuhnBimodal(Bonds,Atoms, options);
        else
            error('Error: distribution type: %s not recognized', dist_type);
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
        
        % Rebuild neighbor lists (choose max_nbr=6 for triangular, or 4 to
        % match your older 4-neighbor layout)
        max_nbr = 6;
        Atoms    = NetworkBuildNeighborLists(Atoms, Bonds, max_nbr);
        
        if strcmp(dist_type,'polydisperse')
            [Nvec] = NetworkGenAssignKuhnPolydisperse(Bonds,Atoms, options);
        elseif strcmp(dist_type,'bimodal')
            [Nvec] = NetworkGenAssignKuhnBimodal(Bonds,Atoms, options);
        else
            error('Error: distribution type: %s not recognized', dist_type);
        end

    else
        error('Unknown network_geometry: %s', options.network_geometry);
    end

    %2. Contruct local density potential
    [LDpot] = NetworkGenConstructLDPotential(Domain, Atoms, Bonds, Nvec, options);
    
    %% D. Scale domain if needed
    if (scale ~= 1.0)
        [Domain, Atoms, Bonds] = NetworkScaleDomain(Domain, Atoms, Bonds, scale);
    end

    %% E. Show visualization and statistics
    NetworkGenVisualize(Domain, Atoms, Bonds, Nvec, scale, options);
    
    %% F. Compute order parameters
    [order] = NetworkComputeOrder(Atoms,Bonds);
    
    %% G. Write data files
    NetworkGenWriteDataFiles(Domain, Atoms, Bonds, Nvec, LDpot, options,order);
    
end

fprintf('Done!\n')
