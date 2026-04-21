try
    addpath('NetworkGen');
    net = network();
    
    % Use small domain for speed
    net.domain.Lx = 20;
    net.domain.Ly = 20;
    net.domain.scale = 1;
    net.domain.write_location = 'tmp_multitype_rulecheck_rerun';
    net.domain.lammps_data_file = 'rulecheck_rerun';
    net.domain.seed = 42;
    
    net.architecture.geometry = 'random';
    net.architecture.rho_atom = 0.5; % dense enough to have many bonds
    
    net.architecture.strand_typology.mode = 'mono';
    net.architecture.strand_typology.mono.n_kuhn = 5;
    
    net.architecture.types.enabled = true;
    net.architecture.types.natom_type = 3;
    net.architecture.types.nbond_type = 2;
    net.architecture.types.atype_mode = 'frac';
    net.architecture.types.btype_mode = 'frac';
    net.architecture.types.atom_frac = [0.5 0.3 0.2];
    net.architecture.types.bond_frac = [0.6 0.4];
    % connectivity: [atomTypeA atomTypeB bondType allowed]
    % 1 2 1 0: Bond Type 1 forbidden between Atom Type 1 and 2
    % 2 3 2 0: Bond Type 2 forbidden between Atom Type 2 and 3
    net.architecture.types.connectivity = [ ...
        1 2 1 0; ...
        1 2 2 1; ...
        2 3 1 1; ...
        2 3 2 0  ...
    ];
    
    net.flags.iplot = false;
    net.flags.ipotential = false;
    net.flags.idefect = false;
    net.flags.isave = true;
    
    net.generateNetwork();
catch ME
    fprintf('Error: %s\n', ME.message);
    % Print stack trace
    for k = 1:length(ME.stack)
        disp(ME.stack(k));
    end
    exit(1);
end
exit(0);
