classdef architecture

properties

    %%% Architecture
    geometry = 'random'                     % 'random', 'hex_lattice', 'bottle_brush'
    strand_typology = assignmentmode();
    types = struct(...
        'enabled',              false, ...
        'natom_type',           1, ...
        'nbond_type',           1, ...
        'atype_mode',           'frac', ...
        'btype_mode',           'frac', ...
        'atom_count',           0, ...
        'bond_count',           0, ...
        'atom_frac',            1, ...
        'bond_frac',            1, ...
        'connectivity',         [ ], ...
        'atype_sel_method',     'random', ...
        'btype_sel_method',     'random', ...
        'btype_same_mol',       2, ...    % bond type for endpoints sharing a molecule ID (btype_sel_method='by_molid')
        'btype_diff_mol',       1 ...     % bond type for endpoints with different molecule IDs
    );
    lattice_spacing =           6;
    spacing_multiplier_mode =   'auto';
    spacing_multiplier =        1.2;
    lattice_disorder_level =     1;
    lattice_disorder_maxfrac =  0.4;
    lattice_max_del_per_node =  1;
    lattice_min_degree_keep =   5;

    rho_atom =                  [];

    %%% Bottle-brush architecture
    %
    % Rod-length dispersity:
    %   rod_dispersity = false -> every rod has exactly Nr atoms.
    %   rod_dispersity = true  -> each rod's atom count is sampled from the
    %                             rod_length assignmentmode. Set its .mode
    %                             ('mono' | 'uniform' | 'poly' | 'bimodal')
    %                             and the matching sub-struct of settings,
    %                             exactly as for strand_typology. The .auto
    %                             flag of rod_length is not used.
    %   Supported rod_length methods: poly -> 'pmf' or 'mono';
    %                                 bimodal -> 'gaussian' or 'single'.
    %   (poly 'geom'/'range' need bond lengths and do not apply to rods.)
    bottlebrush = struct(...
        'Nr',                       10, ...    % baseline atoms per rod (used when rod_dispersity = false)
        'sigma_r',                  1.0, ...   % rod atom separation distance
        'sigma_c_rod',              1.0, ...   % rod atom size (for LD potential)
        'rod_collision_multiplier', 2.0, ...   % rod center separation = multiplier * sigma_r
        'rod_dispersity',           false, ... % when true, per-rod atom count is sampled from rod_length
        'rod_length',               assignmentmode() ...  % distribution governing per-rod atom count
    );

end

end
