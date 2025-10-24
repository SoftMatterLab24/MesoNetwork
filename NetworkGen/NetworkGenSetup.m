function [Domain,Settings] = NetworkGenSetup(options,advancedOptions)

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
node_scatter_max_tries          = Max_atom;      % global limit for scatter loop was 10 * Max_atom
max_tries_per_node_sample       = 200;           % per-node rejection limit

% Bond creation guards
bond_global_try_limit          = 200 * Max_bond; % absolute ceiling
max_attempts_without_progress  = 10  * Max_atom; % local stall guard

% Pruning rule
min_degree_keep    = 1;        % delete nodes with degree <= 1

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
%% B. Sytle specific Parameters
if strcmp(dist_type,'polydisperse')
    % Polydisperse parameters


    if (advancedOptions.iadvancedoptions)

    else
        % use default values

    end

    % Collect params into Domain structure


elseif strcmp(dist_type,'bimodal')
    % Bimodal parameters
    
    N1 = options.bimodal.N1;
    N2 = options.bimodal.N2;

    r1_avg = b*sqrt(N1);
    r2_avg = b*sqrt(N2);

    if (advancedOptions.iadvancedoptions)

    else 
        % use default values
        % Calculate width of inital bond pre-stretch (end-to-end distance)
        sig1 = 0.4*r1_avg;      % standard deviation of bond length distribution
        FWHM1 = 2.355*sig1;     % full width at half maximum 

        sig2 = 1.5*sig1;
        FWHM2 = 2.355*sig2;

        P2 = 0.2;               %Probability constraint to reduce number of long bonds
    end

    % Collect params into Domain structure

end

% Node scattering (min distance)
min_node_sep               = 6*b;            % HARD minimum spacing between any 2 nodes 0.30
min_node_sep2              = min_node_sep^2;    % compare squared dists

Domain.min_node_sep        = min_node_sep;      
Domain.min_node_sep2       = min_node_sep2;
end