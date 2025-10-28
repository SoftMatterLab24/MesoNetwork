function NetworkGenWriteDataFiles(Domain, Atoms, Bonds, options,Nvec)
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

% ---------- Prep paths ----------
if isfield(options,'write_location') && ~isempty(options.write_location)
    outdir = options.write_location;
else
    outdir = '.';
end
if ~exist(outdir,'dir'), mkdir(outdir); end

if ~isfield(options,'lammps_data_file') || isempty(options.lammps_data_file)
    error('options.lammps_data_file is required.');
end
data_path = fullfile(outdir, options.lammps_data_file);

if ~isfield(options,'bond_table_file') || isempty(options.bond_table_file)
    bond_table_file = 'bond.table';
else
    bond_table_file = options.bond_table_file;
end
bondtable_path = fullfile(outdir, bond_table_file);

% ---------- Gather counts & domain ----------
Atom_count = size(Atoms,1);
Bond_count = size(Bonds,1);

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
fprintf(fid, '%d atom types\n', 1);
fprintf(fid, '%d bond types\n', 1);
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
    fprintf(fid, '%d 1 %d %d\n', Bonds(i,1), Bonds(i,2), Bonds(i,3));
end

fclose(fid);
fprintf('Wrote %s with %d atoms and %d bonds.\n', data_path, Atom_count, Bond_count);

% ---------- Write bond.table ----------
b = options.b;  % Kuhn length
if Bond_count > 0
    % Use provided Nvec if available; otherwise compute a fallback.
    if isempty(Nvec)
        if isfield(options,'dist_type') && strcmpi(options.dist_type,'polydisperse')
            Nvec = NetworkGenAssignKuhnPolydisperse(Bonds, options);
        else
            % TODO: plug a bimodal assigner later; for now default to ones.
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
    fprintf('Wrote %s with %d entries. b = %.6g\n', bondtable_path, Bond_count, b);
else
    % Still create an empty bond.table with header for consistency
    fidBT = fopen(bondtable_path,'w');
    if fidBT >= 0
        fprintf(fidBT, '# Chain stats\n\n');
        fprintf(fidBT, 'KEY\n');
        fprintf(fidBT, 'N 0\n\n');
        fclose(fidBT);
    end
    fprintf('No bonds; wrote empty %s.\n', bondtable_path);
end
