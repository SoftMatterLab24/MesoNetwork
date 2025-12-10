function [Domain] = NetworkGenSetup(options,advancedOptions)
% NetworkGenSetup - Setup the simulation domain and parameters
%
% INPUT:
%   options: structure with general options
%   advancedOptions: structure with advanced options
% OUTPUT:
%   Domain: structure with domain and parameter settings

%% Extract relevant options
Lx = options.Lx;
Ly = options.Ly;
b  = options.b;

dist_type = options.dist_type;

%% A. General Domain Parameters
% Domain size 
xlo = -Lx*b; xhi = Lx*b;
ylo = -Ly*b; yhi = Ly*b;
zlo = -20*b; zhi = 20*b;

% Network Size Caps
Rho_atom           = 0.0078;                                   % intial atom density
Max_atom           = ceil(Rho_atom*(xhi-xlo)*(yhi-ylo));       % requested # nodes
Max_peratom_bond   = 5;                                        % degree cap per node (fits neighbor slots)
Max_bond           = round(0.5*Max_atom*Max_peratom_bond);     % cap on bonds

% Atom creation guards
node_scatter_max_tries          = Max_atom;      % global limit for scatter loop
max_tries_per_node_sample       = 200;           % per-node rejection limit

% Bond creation guards
bond_global_try_limit          = 200 * Max_bond; % absolute ceiling
max_attempts_without_progress  = 10  * Max_atom; % local stall guard

% Pruning rule
min_degree_keep    = 1;        % delete nodes with degree <= 1

% Node scattering (min distance)
min_node_sep               = 6*b;               % HARD minimum spacing between any 2 nodes 0.30
min_node_sep2              = min_node_sep^2;    % compare squared dists

% Write to Domain structure
Domain.xlo = xlo; Domain.xhi = xhi;
Domain.ylo = ylo; Domain.yhi = yhi;
Domain.zlo = zlo; Domain.zhi = zhi;

Domain.Max_atom                 = Max_atom;
Domain.Max_bond                 = Max_bond;
Domain.Max_peratom_bond         = Max_peratom_bond;

Domain.node_scatter_max_tries   = node_scatter_max_tries;
Domain.bond_global_try_limit    = bond_global_try_limit;
Domain.max_attempts_without_progress = max_attempts_without_progress;

Domain.max_tries_per_node_sample  = max_tries_per_node_sample;
Domain.min_degree_keep            = min_degree_keep;

Domain.min_node_sep        = min_node_sep;      
Domain.min_node_sep2       = min_node_sep2;

%% B. Override advanced options
if (advancedOptions.iadvancedoptions)

   Rho_atom = advancedOptions.Rho_atom;
   Max_peratom_bond = advancedOptions.Max_peratom_bond;
   bond_global_try_limit_multiplier = advancedOptions.bond_global_try_limit_multiplier;
   max_attempts_without_progress_multiplier = advancedOptions.max_attempts_without_progress_multiplier;
   min_degree_keep = advancedOptions.min_degree_keep;
   cutoff_multiply = advancedOptions.cutoff_multiply; %units of b

   % Recalculate
    % Network Size Caps
    Max_atom           = ceil(Rho_atom*(xhi-xlo)*(yhi-ylo));       % requested # nodes
    Max_bond           = round(0.5*Max_atom*Max_peratom_bond);     % cap on bonds

    % Atom creation guards
    node_scatter_max_tries          = Max_atom;      % global limit for scatter loop
    max_tries_per_node_sample       = 200;           % per-node rejection limit

    % Bond creation guards
    bond_global_try_limit          = bond_global_try_limit_multiplier * Max_bond; % absolute ceiling
    max_attempts_without_progress  = max_attempts_without_progress_multiplier  * Max_atom; % local stall guard

    % Node scattering (min distance)
    min_node_sep               = cutoff_multiply*b; % HARD minimum spacing between any 2 nodes 0.30
    min_node_sep2              = min_node_sep^2;    % compare squared dists

    % Recollect parameters into Domain struct
    Domain.xlo = xlo; Domain.xhi = xhi;
    Domain.ylo = ylo; Domain.yhi = yhi;
    Domain.zlo = zlo; Domain.zhi = zhi;

    Domain.Max_atom                 = Max_atom;
    Domain.Max_bond                 = Max_bond;
    Domain.Max_peratom_bond         = Max_peratom_bond;

    Domain.node_scatter_max_tries   = node_scatter_max_tries;
    Domain.bond_global_try_limit    = bond_global_try_limit;
    Domain.max_attempts_without_progress = max_attempts_without_progress;

    Domain.max_tries_per_node_sample  = max_tries_per_node_sample;
    Domain.min_degree_keep            = min_degree_keep;

    Domain.min_node_sep        = min_node_sep;      
    Domain.min_node_sep2       = min_node_sep2;
else
    % continue to use defaults
end

%% C. Update file names
data_file = options.lammps_data_file;

type_prefix = '';

end