function [] = SetupDomain(obj)
% -------------------------------------------------------------------------
% SetupDomain
% - Builds derived domain/runtime quantities from the network object
% - Reads user settings from obj.domain / obj.architecture / obj.peratom
% - Writes derived values back into obj.domain
%
% INPUT:
%   obj : network class object
%
% OUTPUT:
%   none
%   (obj.domain is modified in place)
% -------------------------------------------------------------------------

    %% ---------------------- Extract primary inputs ----------------------
    Lx = obj.domain.Lx;
    Ly = obj.domain.Ly;
    Lz = obj.domain.Lz;
    b  = obj.domain.b;

    %% ---------------------- General domain bounds -----------------------
    obj.domain.xlo = -Lx*b;
    obj.domain.xhi =  Lx*b;

    obj.domain.ylo = -Ly*b;
    obj.domain.yhi =  Ly*b;

    obj.domain.zlo = -Lz*b;
    obj.domain.zhi =  Lz*b;

    %% ---------------------- Default / user controls ---------------------
    % Per-atom bond cap
    if isempty(obj.peratom.Max_peratom_bond)
        Max_peratom_bond = 5;
    else
        Max_peratom_bond = obj.peratom.Max_peratom_bond;
    end

    % Pruning rule
    if isempty(obj.peratom.min_degree_keep)
        min_degree_keep = 1;
    else
        min_degree_keep = obj.peratom.min_degree_keep;
    end

    % Node spacing: use architecture.lattice_spacing directly
    if isempty(obj.architecture.lattice_spacing)
        min_node_sep = 6*b;
    else
        min_node_sep = obj.architecture.lattice_spacing * b;
    end
    min_node_sep2 = min_node_sep^2;

    % Atom density
    if ~isempty(obj.architecture.rho_atom)
        Rho_atom = obj.architecture.rho_atom;
    else
        Rho_atom = 0.7/((pi/4)*min_node_sep2);
        obj.architecture.rho_atom = Rho_atom;
    end
    
    %% ---------------------- Network size caps ---------------------------
    obj.domain.Max_atom = ceil( ...
        Rho_atom * ...
        (obj.domain.xhi - obj.domain.xlo) * ...
        (obj.domain.yhi - obj.domain.ylo) );

    obj.peratom.Max_peratom_bond = Max_peratom_bond;

    obj.domain.Max_bond = round( ...
        0.5 * obj.domain.Max_atom * obj.peratom.Max_peratom_bond );

    %% ---------------------- Atom creation guards ------------------------
    obj.domain.node_scatter_max_tries    = obj.domain.Max_atom;
    obj.domain.max_tries_per_node_sample = 200;

    %% ---------------------- Bond creation guards ------------------------
    obj.domain.bond_global_try_limit         = 200 * obj.domain.Max_bond;
    obj.domain.max_attempts_without_progress = 10  * obj.domain.Max_atom;

    %% ---------------------- Store spacing / pruning ---------------------
    obj.domain.min_node_sep    = min_node_sep;
    obj.domain.min_node_sep2   = min_node_sep2;



end