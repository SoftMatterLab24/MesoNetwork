%% Create a network object
n = network;

n.dist_type = 'bimodal';

% Number of networks to generate
n.Nreplicates = 1;

% Kuhn length
n.b = 1.6;     % Kuhn length (in nm)

% Domain size
n.Lx = 100;    % Domain size in x (in units of b)
n.Ly = 100;    % Domain size in y (in units of b)

% Seed options
n.imanualseed = false;  % true: manual seed; false: random seed
n.seed = [1];

% Visualization
% n.iplot = true;    % Show 

% Save options
% n.isave = true;  % Save data files
% n.lammps_data_file   = 'PronyNetwork_nano_1300x800.dat';          % Prefix file name for LAMMPS data output
% n.lammps_visual_file = 'PronyVisual_10000_nano_1300x800.dat';     % Prefix file name for LAMMPS visualization output
% n.bond_table_file    = 'bond.table';                              % File name for bond table output   
% n.write_location     = './networks';                              % Location to write output files

%% --------------------- Polydisperse Options ----------------------

n.distribution_assignment_mode_poly = 'pmf';  % Kuhn segment assigment method: 'geom' | 'range' | 'pmf'

%% --------------------- Bimodal Options ---------------------------
n.N1 = 20; 
n.N2 = 250;

n.bin_window_method = 'manual';             % Method for determining bin width of bimodal dist: 'manual' or 'adaptive'               

n.distribution_assignment_mode = 'gaussian';  % Kuhn segment assigment method: 'single' or 'geom' or 'gaussian'
n.distribution_height_mode = 'fixed';         % Distribution height method: 'prob' or 'fixed'
n.long_first = true;                          % enable long-first mode

% Height mode settings (only one is used)
n.P = 0.2;        % Prob: desired fraction of type 2 bonds
n.N2_bonds = 0; % Fixed: desired number of type 2 bonds

% Manual Window parameters (if bin_window_method = manual)
n.dr1 = 5;
n.dr2 = 5;

n.iadvancedoptions = false;

%% Make the network
[Domain,Atoms,Bonds,Nvec] = generateNetwork(n);

 % extract data from atoms
            x = Atoms(:,2);
            y = Atoms(:,3);
            Total_bond = length(Bonds);
            
            figure; hold on
            scatter(x,y, 6, 'k', 'filled')

            for k = 1:Total_bond
                k;
                i1 = Bonds(k,2);
                i2 = Bonds(k,3);
                
                if Bonds(k,5) == 1
                    plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'k-');
                else
                    plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'r-');
                end
            end