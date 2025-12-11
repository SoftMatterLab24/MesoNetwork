function NetworkGenWriteDataFiles(Domain,Atoms,Bonds,Nvec,LDpot,options)
%NetworkGenWriteDataFiles - Write LAMMPS data file and bond.table
% 
% Files:
%   - <write_location>/<lammps_data_file>
%   - <write_location>/<bond_table_file>
%
% LAMMPS atoms layout written as: 
%   atomID  molID  atomType  diameter  density   x  y  z
% Bonds layout written as:
%   bondID  bondType  atom1  atom2
%
% bond.table layout:
%   # Chain stats
%   KEY
%   N <#bonds>
%   <bondID> <id1> <id2> <N> <b>

if options.isave ~= true
    warning('Did not write data files because options.isave = false.');
    return; % skip writing files
end

fprintf('   Writing network to data file...\n');

% ---------- Prep filenames ---------- 

% Determine network type
if strcmpi(options.dist_type,'bimodal')
    networkTY_prefix = 'BD';
elseif strcmpi(options.dist_type,'polydisperse') && ~strcmpi(options.polydisperse.distribution_assignment_mode,'mono')
    networkTY_prefix = 'PD';
else
    networkTY_prefix = 'MD';
end

if options.save_name_mode == true
    % Auto naming mode: add sample info to file names
    sample_label = sprintf('%s_%s_%s', ...
        networkTY_prefix, ...
        options.sample_suffix, ...
        options.replicate_suffix);
       
    options.lammps_data_file   = sprintf('%s_%s.dat',options.lammps_data_file, sample_label);
    options.lammps_visual_file = sprintf('%s_%s.dat',options.lammps_visual_file, sample_label);
    options.bond_table_file    = sprintf('%s_%s.table',options.bond_table_file, sample_label);
    options.log_file           = sprintf('%s.log', sample_label);
    options.pot_file           = sprintf('%s.LD.table',sample_label);
else
    % Fixed names mode: use provided names directly
    options.lammps_data_file   = sprintf('%s.dat',options.lammps_data_file);
    options.lammps_visual_file = sprintf('%s.dat',options.lammps_visual_file);
    options.bond_table_file    = sprintf('%s.table',options.bond_table_file);
    options.log_file           = sprintf('%s.log',options.lammps_data_file);
    options.pot_file           = sprintf('%s.LD.table',options.lammps_data_file);
end

% ---------- Prep paths ----------
% Output directory
if isfield(options,'write_location') && ~isempty(options.write_location)
    outdir = options.write_location;
else
    outdir = '.';
end
if ~exist(outdir,'dir'), mkdir(outdir); end

% LAMMPS data file path
if ~isfield(options,'lammps_data_file') || isempty(options.lammps_data_file)
    error('options.lammps_data_file is required.');
end
data_path = fullfile(outdir, options.lammps_data_file);

% bond.table path
if ~isfield(options,'bond_table_file') || isempty(options.bond_table_file)
    bond_table_file = 'bond.table';
else
    bond_table_file = options.bond_table_file;
end
bondtable_path = fullfile(outdir, bond_table_file);

% manybody potential path
potfile_path = fullfile(outdir, options.pot_file);

% log file path
logfile_path = fullfile(outdir,options.log_file);

% ---------- Gather counts & domain ----------
Atom_count = size(Atoms,1);
Bond_count = size(Bonds,1);

natype = 1;

[~, ncol] = size(Bonds);
if ncol > 4
    nbtype = max(Bonds(:,5));
else
    nbtype = 1;
end

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;
zlo = Domain.zlo; zhi = Domain.zhi;

% ---------- Write LAMMPS data file ----------
fid = fopen(data_path,'w');
if fid < 0
    error('Could not open %s for writing.', data_path);
end

fprintf(fid, '\n\n');
fprintf(fid, '%d atoms\n', Atom_count);
fprintf(fid, '%d bonds\n', Bond_count);
fprintf(fid, '%d atom types\n', natype);
fprintf(fid, '%d bond types\n', nbtype);
fprintf(fid, '%.16g %.16g xlo xhi\n', xlo, xhi);
fprintf(fid, '%.16g %.16g ylo yhi\n', ylo, yhi);
fprintf(fid, '%.16g %.16g zlo zhi\n', zlo, zhi);
fprintf(fid, '\n');

fprintf(fid, 'Atoms #bpm/sphere\n\n');
% Format: atomID  molID  atomType  diameter  density   x  y  z
for i = 1:Atom_count
    fprintf(fid, '%d 1 1 1 1 %.16g %.16g %.16g\n', ...
        Atoms(i,1), Atoms(i,2), Atoms(i,3), Atoms(i,4));
end

fprintf(fid, '\nBonds\n\n');
% Format: bondID  bondType  atom1  atom2
for i = 1:Bond_count
    if ncol > 4; btype = Bonds(i,5); else btype = 1; end
        
    fprintf(fid, '%d %d %d %d\n', Bonds(i,1), btype, Bonds(i,2), Bonds(i,3));
end

fclose(fid);
fprintf('   Wrote %s with %d atoms and %d bonds.\n', data_path, Atom_count, Bond_count);

% ---------- Write bond.table ----------
b = options.b;  % Kuhn length
if Bond_count > 0
    % Use provided Nvec if available; otherwise compute a fallback.
    if isempty(Nvec)
        if isfield(options,'dist_type') && strcmpi(options.dist_type,'polydisperse')
            Nvec = NetworkGenAssignKuhnPolydisperse(Bonds, options);
        elseif isfield(options,'dist_type') && strcmpi(options.dist_type,'bimodal')
            Nvec = NetworkGenAssignKuhnBimodal(Bonds, options);
        else
            % fallback
            Nvec = ones(Bond_count,1);
        end
    end

    % Safety: shape check
    if ~isvector(Nvec) || numel(Nvec) ~= Bond_count
        error('Nvec must be a vector of length Bond_count (= %d).', Bond_count);
    end
    Nvec = Nvec(:);  % column

    fidBT = fopen(bondtable_path,'w');
    if fidBT < 0
        error('Could not open %s for writing.', bondtable_path);
    end

    fprintf(fidBT, '# Chain stats\n\n');
    fprintf(fidBT, 'KEY\n');
    fprintf(fidBT, 'N %d\n\n', Bond_count);

    % Lines: bond-id   i1-atom   i2-atom   N    b
    for k = 1:Bond_count
        fprintf(fidBT, '%d %d %d %d %.8g\n', ...
            Bonds(k,1), Bonds(k,2), Bonds(k,3), Nvec(k), b);
    end

    fclose(fidBT);
    fprintf('   Wrote %s with %d entries. b = %.6g\n', bondtable_path, Bond_count, b);
else
    % Still create an empty bond.table with header for consistency
    fidBT = fopen(bondtable_path,'w');
    if fidBT >= 0
        fprintf(fidBT, '# Chain stats\n\n');
        fprintf(fidBT, 'KEY\n');
        fprintf(fidBT, 'N 0\n\n');
        fclose(fidBT);
    end
    fprintf('   No bonds; wrote empty %s.\n', bondtable_path);


    % Print out the appropriate LJ lengthscale for this specific network
    Total_kuhn_segment = sum(Nvec);
    Node_assigned_radius =(Total_kuhn_segment/Atom_count)^(1/2);
    sigma = Node_assigned_radius / (2^(1/6));                      % convert to LJ sigma parameter
    fprintf('The 6-12 Lennard-Jones radius is equal to %g\n',sigma);
end

% LD pot is struct with fields:
% LDpot.N_LD          = 1;        % single local density potential
% LDpot.N_rho         = N_rho;    % number of density points
% LDpot.R_lower       = R1;       % 10% overlap with other repulisve potentials 
% LDpot.R_upper       = R2;       % outer radius of density calculation
% LDpot.rc            = rc;       % cutoff radius for bpm/spring repulsion
% LDpot.rho_min       = rho_min; 
% LDpot.rho0          = rho0;
% LDpot.rho_max       = rho_max;
% LDpot.drho          = drho;
% LDpot.pot_density   = pot_density;

% ---------- Write manybody pot file ----------
fidP = fopen(potfile_path,'w');

if fidP < 0 
    error('Could not open %s for writing.',potfile_path);
end

fprintf(fidP, '\n\n');
fprintf(fidP, '%d %d \n',LDpot.N_LD,LDpot.N_rho);
fprintf(fidP, '\n');
fprintf(fidP, '%2.4f %2.4f \n',LDpot.R_lower,LDpot.R_upper);
fprintf(fidP, '1 \n');
fprintf(fidP, '1 \n');
fprintf(fidP, '%2.4f %2.4f %2.4f \n',LDpot.rho_min,LDpot.rho_max,LDpot.drho);

for i = 1:length(LDpot.pot_density)
    fprintf(fidP, '%2.8f \n', LDpot.pot_density(i));
end

fclose(fidP);


% ---------- Write log file ----------
fidL = fopen(logfile_path,'w');
if fidL < 0
    error('Could not open %s for writing.','network_log.txt');
end

% Need the min and max 
xAlo = min(Atoms(:,2));
xAhi = max(Atoms(:,2));
yAlo = min(Atoms(:,3));
yAhi = max(Atoms(:,3));

fprintf(fidL,"Sample type:      %s\n", networkTY_prefix);
fprintf(fidL,"Sample number:    %d\n", sscanf(options.sample_suffix,'SMP%d'));
fprintf(fidL,"Replicate number: %d\n", sscanf(options.replicate_suffix,'N%d'));
fprintf(fidL,"Number of atoms:  %d\n", Atom_count);
fprintf(fidL,"Number of bonds:  %d\n", Bond_count);
fprintf(fidL,"Domain size: \n xlo xhi: %.4f %.4f \n ylo yhi: %.4f %.4f \n zlo zhi: %.4f %.4f\n", xAlo, xAhi, yAlo, yAhi, zlo, zhi);
fprintf(fidL,"Equilibrium density:           %.6f\n", LDpot.rho0);
fprintf(fidL,"Lower cutoff radius (R_lower): %.4f\n", LDpot.R_lower/options.b);
fprintf(fidL,"Upper cutoff radius (R_upper): %.4f\n", LDpot.R_upper/options.b);
fprintf(fidL,"BPM cutoff radius (rc):        %.4f\n", LDpot.rc/options.b);
fprintf(fidL,"Crosslink density:             %.6f\n", Atom_count / ((xhi - xlo)*(yhi - ylo)*(zhi - zlo)));
fprintf(fidL,"Kuhn segments per crosslink:   %.4f\n", sum(Nvec) / Atom_count);

fclose(fidL);
fprintf('   Wrote %s\n',logfile_path)