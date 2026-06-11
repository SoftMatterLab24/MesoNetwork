function WriteDataFiles(obj, Atoms, Bonds, Nvec, LDpot, TypeData)
% -------------------------------------------------------------------------
% WriteDataFiles
% - Write LAMMPS data file
% - Write bond table
% - Optionally write local-density potential table
%
% Atoms column layout:
%   [ID | molID | X | Y | Z | deg | nbrs...]
%   molID is always present (set by every AddAtoms* routine).
%
% INPUTS
%   obj      : network object
%   Atoms    : atom array
%   Bonds    : bond array [bondID | id1 | id2 | L0 | type]
%   Nvec     : per-bond Kuhn segment counts
%   LDpot    : local density potential struct, or []
%   TypeData : runtime type-label struct from AssignMultiType, or []
%              (atom types). Bond types are already written in Bonds(:,5).
% -------------------------------------------------------------------------

    if nargin < 6
        TypeData = [];
    end

    if ~obj.flags.isave
        obj.log.print('   Did not write data files because obj.flags.isave = false.\n');
        return;
    end

    obj.log.print('   Writing network data files...\n');

    outdir = obj.domain.write_location;
    if isempty(outdir)
        outdir = '.';
    end

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    data_path      = fullfile(outdir, obj.log.lammps_data_file);
    bondtable_path = fullfile(outdir, obj.log.bond_table_file);

    if ~isempty(obj.log.pot_file)
        potfile_path = fullfile(outdir, obj.log.pot_file);
    else
        potfile_path = '';
    end

    Atom_count = size(Atoms,1);
    Bond_count = size(Bonds,1);

    xlo = obj.domain.xlo; xhi = obj.domain.xhi;
    ylo = obj.domain.ylo; yhi = obj.domain.yhi;
    zlo = obj.domain.zlo; zhi = obj.domain.zhi;

    % ---------------------------------------------------------------------
    % Molecule IDs: always at col 2 of Atoms
    % ---------------------------------------------------------------------
    if size(Atoms, 2) >= 2
        mol_id_vec = Atoms(:, 2);
        mol_id_vec(mol_id_vec == 0) = 1;
    else
        mol_id_vec = ones(Atom_count, 1);
    end

    if numel(mol_id_vec) ~= Atom_count
        error('WriteDataFiles: molecule ID vector must have length %d.', Atom_count);
    end

    % ---------------------------------------------------------------------
    % Atom / bond type vectors
    % Bond types come directly from Bonds(:,5), which AssignMultiType has
    % populated by this point (with a safe default of 1 when typing is off).
    % ---------------------------------------------------------------------
    atom_type_vec = ones(Atom_count, 1);

    if ~isempty(TypeData) && isstruct(TypeData) && isfield(TypeData, 'enabled') && TypeData.enabled
        if isfield(TypeData, 'atom_types') && ~isempty(TypeData.atom_types)
            atom_type_vec = TypeData.atom_types(:);
        end
        natype = max(max(atom_type_vec), TypeData.natom_type);
    else
        natype = 1;
    end

    if size(Bonds,2) >= 5 && ~isempty(Bonds)
        bond_type_vec = Bonds(:,5);
    else
        bond_type_vec = ones(Bond_count, 1);
    end

    if ~isempty(TypeData) && isstruct(TypeData) && isfield(TypeData, 'enabled') && TypeData.enabled
        if isempty(bond_type_vec)
            nbtype = max(1, TypeData.nbond_type);
        else
            nbtype = max(max(bond_type_vec), TypeData.nbond_type);
        end
    else
        if isempty(bond_type_vec)
            nbtype = 1;
        else
            nbtype = max(max(bond_type_vec), 1);
        end
    end

    if numel(atom_type_vec) ~= Atom_count
        error('WriteDataFiles: atom type vector must have length %d.', Atom_count);
    end

    if numel(bond_type_vec) ~= Bond_count
        error('WriteDataFiles: bond type vector must have length %d.', Bond_count);
    end

    % ---------------------------------------------------------------------
    % Write LAMMPS data file
    % ---------------------------------------------------------------------
    fid = fopen(data_path, 'w');
    if fid < 0
        error('WriteDataFiles: cannot open %s for writing.', data_path);
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
    % atomID molID atomType diameter density x y z
    %   coords at cols 3 (x), 4 (y), 5 (z) under new layout
    for i = 1:Atom_count
        fprintf(fid, '%d %d %d 1 1 %.16g %.16g %.16g\n', ...
            Atoms(i,1), mol_id_vec(i), atom_type_vec(i), ...
            Atoms(i,3), Atoms(i,4), Atoms(i,5));
    end

    fprintf(fid, '\nBonds\n\n');
    for i = 1:Bond_count
        btype = bond_type_vec(i);
        fprintf(fid, '%d %d %d %d\n', Bonds(i,1), btype, Bonds(i,2), Bonds(i,3));
    end

    fclose(fid);

    obj.log.print('   Wrote %s with %d atoms and %d bonds.\n', ...
        data_path, Atom_count, Bond_count);

    % ---------------------------------------------------------------------
    % Write bond.table
    % ---------------------------------------------------------------------
    b = obj.domain.b;

    if isempty(Nvec)
        Nvec = ones(Bond_count,1);
    end

    if ~isvector(Nvec) || numel(Nvec) ~= Bond_count
        error('WriteDataFiles: Nvec must be a vector of length %d.', Bond_count);
    end
    Nvec = Nvec(:);

    fidBT = fopen(bondtable_path, 'w');
    if fidBT < 0
        error('WriteDataFiles: cannot open %s for writing.', bondtable_path);
    end

    fprintf(fidBT, '# Chain stats\n\n');
    fprintf(fidBT, 'KEY\n');
    fprintf(fidBT, 'N %d\n\n', Bond_count);

    for k = 1:Bond_count
        fprintf(fidBT, '%d %d %d %d %.8g\n', ...
            Bonds(k,1), Bonds(k,2), Bonds(k,3), Nvec(k), b);
    end

    fclose(fidBT);

    obj.log.print('   Wrote %s with %d entries. b = %.6g\n', ...
        bondtable_path, Bond_count, b);

    % ---------------------------------------------------------------------
    % Optional manybody potential file
    % ---------------------------------------------------------------------
    if obj.flags.ipotential && ~isempty(LDpot)

        fidP = fopen(potfile_path, 'w');
        if fidP < 0
            error('WriteDataFiles: cannot open %s for writing.', potfile_path);
        end

        fprintf(fidP, '\n\n');
        fprintf(fidP, '%d %d \n', LDpot.N_LD, LDpot.N_rho);
        fprintf(fidP, '\n');
        fprintf(fidP, '%2.4f %2.4f \n', LDpot.R_lower, LDpot.R_upper);
        fprintf(fidP, '1 \n');
        fprintf(fidP, '1 \n');
        fprintf(fidP, '%2.4f %2.4f %2.8f \n', ...
            LDpot.rho_min, LDpot.rho_max, LDpot.drho);

        for i = 1:length(LDpot.pot_density)
            fprintf(fidP, '%2.8f \n', LDpot.pot_density(i));
        end

        fclose(fidP);

        if isfield(LDpot, 'type')
            obj.log.print('   Wrote %s (%s LD potential)\n', potfile_path, LDpot.type);
        else
            obj.log.print('   Wrote %s\n', potfile_path);
        end
    elseif obj.flags.ipotential
        obj.log.print('   Skipped potential write because LDpot is empty.\n');
    end

end
