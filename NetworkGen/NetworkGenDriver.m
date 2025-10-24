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

%% --------------------------- Global settings ----------------------
% Distribution type: 'bimodal' or 'polydisperse'
dist_type = 'bimodal';

% Number of networks to generate
Nreplicates = 1;

% Kuhn length
b = 1.6;     % Kuhn length (in nm)

% Domain size
Lx = 800;    % Domain size in x (in units of b)
Ly = 800;    % Domain size in y (in units of b)

% Seed options
imanualseed = true;  % true: manual seed; false: random seed
seed = [1];

% Visualization
iplot = true;    % Show 

% Save options
lammps_data_file   = 'PronyNetwork_nano_1300x800.dat';          % Prefix file name for LAMMPS data output
lammps_visual_file = 'PronyVisual_10000_nano_1300x800.dat';     % Prefix file name for LAMMPS visualization output
bond_table_file    = 'bond.table';                              % File name for bond table output   
write_location     = './networks';                              % Location to write output files

%% --------------------- Polydisperse Options ----------------------



%% --------------------- Bimodal Options ---------------------------
N1 = 1; 
N2 = 1;

P = 0.5;    % desired fraction of type 2 bonds

%% --------------------- Advanced Options ---------------------------
iadvancedoptions = false;

% Rcut = 10*b;    % cutoff distance for bonding (in units of b)

%% --------------------- Network Generation ------------------------
% DO NOT EDIT BELOW THIS LINE

% prepare options structure
options.dist_type          = dist_type;
options.Nreplicates        = Nreplicates;
options.b                  = b;
options.Lx                 = Lx;
options.Ly                 = Ly;
options.imanualseed        = imanualseed;
options.seed               = seed;
options.iplot              = iplot;
options.lammps_data_file   = lammps_data_file;
options.lammps_visual_file = lammps_visual_file;
options.bond_table_file    = bond_table_file;
options.write_location     = write_location;

% Polydisperse options
% (none for now)
% options.polydisperse.<param> = <value>;

% Bimodal options
options.bimodal.N1                 = N1;
options.bimodal.N2                 = N2;
options.bimodal.P                  = P;

% Additional advanced options
advancedOptions.iadvancedoptions = iadvancedoptions;
if iadvancedoptions
    % (add advanced options here)
end

for ii = 1:Nreplicates
    fprintf('Generating network replicate %d of %d...\n',ii,Nreplicates);
    

    %% A. Prepare replicate-specific information
    % 1. Set seed
    if imanualseed
        if length(seed) ~= Nreplicates
            error('Error: not enough manual seeds provided for the number of replicates');
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
        [Bonds] = NetworkGenConnectNodesPolydisperse(Domain,Atoms,options);
    elseif strcmp(dist_type,'bimodal')
        % Bimodal network
        [Bonds] = NetworkGenConnectNodesBimodal(Domain,Atoms,options);
    else
        error('Error: distribution type not recognized');
    end

    %% E. Show visualization and statistics
    NetworkGenVisualize(Domain,Atoms,Bonds,options);

    %% F. Write data files
    NetworkGenWriteDataFiles(Domain,Atoms,Bonds,options);
    
end

