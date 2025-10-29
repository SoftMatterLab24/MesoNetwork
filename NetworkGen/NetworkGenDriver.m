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
% Distribution type: 'bimodal' or 'polydisperse'
dist_type = 'bimodal';

% Number of networks to generate
Nreplicates = 1;

% Kuhn length
b = 1.6;     % Kuhn length (in nm)

% Domain size
Lx = 100;    % Domain size in x (in units of b)
Ly = 100;    % Domain size in y (in units of b)

% Seed options
imanualseed = true;  % true: manual seed; false: random seed
seed = [1];

% Visualization
iplot = true;    % Show 

% Save options
isave = true;  % Save data files
lammps_data_file   = 'PronyNetwork_nano_1300x800.dat';          % Prefix file name for LAMMPS data output
lammps_visual_file = 'PronyVisual_10000_nano_1300x800.dat';     % Prefix file name for LAMMPS visualization output
bond_table_file    = 'bond.table';                              % File name for bond table output   
write_location     = './networks';                              % Location to write output files

%% --------------------- Polydisperse Options ----------------------

distribution_assignment_mode = 'pmf';  % Kuhn segment assigment method: 'geom' | 'range' | 'pmf'

%% --------------------- Bimodal Options ---------------------------
N1 = 50; 
N2 = 250;

distribution_assignment_mode = 'single';    % Kuhn segment assigment method: 'single' or 'geom'
distribution_height_mode = 'fixed';          % Distribution height method: 'prob' or 'fixed'

% Height mode settings (only one is used)
P = 0.2;        % Prob: desired fraction of type 2 bonds
N2_bonds = 500; % Fixed: desired number of type 2 bonds

%% --------------------- Advanced Options --------------------------
iadvancedoptions = false;

%% --------------------- Network Generation ------------------------
% !!!DO NOT EDIT BELOW THIS LINE!!!

% prepare options structure
options.dist_type          = dist_type;
options.Nreplicates        = Nreplicates;
options.b                  = b;
options.Lx                 = Lx;
options.Ly                 = Ly;
options.imanualseed        = imanualseed;
options.seed               = seed;
options.iplot              = iplot;
options.isave              = isave;
options.lammps_data_file   = lammps_data_file;
options.lammps_visual_file = lammps_visual_file;
options.bond_table_file    = bond_table_file;
options.write_location     = write_location;

% A. Polydisperse options
% ------------------------------------------------------------------
% --- mode selection ---
options.polydisperse.distribution_assignment_mode = distribution_assignment_mode;   % 'geom' | 'range' | 'pmf'

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
options.bimodal.N1                 = N1;
options.bimodal.N2                 = N2;

% --- mode selection ---
% 'single' mode: applies N1 and N2 directly based on geometry (N ≈ (L/b)^2)
% 'range' mode: uses: b (global), kuhn_rounding, min_N
options.bimodal.distribution_assignment_mode = distribution_assignment_mode; % 'single' | 'geom'

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

    % 2. Set replicate-specific file names
    if Nreplicates > 1
        options.lammps_data_file   = sprintf('%04d_%s',ii,lammps_data_file);
        options.lammps_visual_file = sprintf('%04d_%s',ii,lammps_visual_file);
        options.bond_table_file    = sprintf('%04d_%s',ii,bond_table_file);
    end

    %% B. Setup the system
    [Domain] = NetworkGenSetup(options,advancedOptions);

    %% C. Scatter nodes with minimum spacing
    [Atoms] = NetworkGenScatterNodes(Domain);

    %% D. Connect nodes within bonds
    if strcmp(dist_type,'polydisperse')
        % Polydisperse network
        [Atoms,Bonds] = NetworkGenConnectNodesPolydisperse(Domain,Atoms,options);
        [Nvec] = NetworkGenAssignKuhnPolydisperse(Bonds, options);
    elseif strcmp(dist_type,'bimodal')
        % Bimodal network
        [Atoms,Bonds] = NetworkGenConnectNodesBimodal(Domain,Atoms,options);
        [Nvec] = NetworkGenAssignKuhnBimodal(Bonds, options);
    else
        error('Error: distribution type: %s not recognized', dist_type);
    end

    %% E. Show visualization and statistics
    NetworkGenVisualize(Domain,Atoms,Bonds,Nvec,options);

    %% F. Write data files
    NetworkGenWriteDataFiles(Domain,Atoms,Bonds,Nvec,options);
    
end

