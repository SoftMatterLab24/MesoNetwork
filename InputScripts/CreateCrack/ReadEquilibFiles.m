clc; clear; close all;

%
%location to lammps outputs (atoms_equilib.dump,bonds_equilib.dump)
loc = 'C:\Users\zwhit\Downloads\polydisperse_net_generator\Runs\bimodal\Unnotched\Run2';
bond_name = 'bonds_equilib.dump';

for ii = 1:5

    bsave = strcat(loc,'\',sprintf('%03d',ii),'\','bimodal_bondequib_',sprintf('%03d',ii),'.dump');
    bfile = strcat(loc,'\',sprintf('%03d',ii),'\',bond_name);
    [bAtom1,bAtom2,bType,bLength,bForce,N,b] = parse_equidump_bond_fun(bfile);

    bondAtom1 = bAtom1{end};
    bondAtom2 = bAtom2{end};
    bondType = bType{end};
    bondLength = bLength{end};
    bondForce = bForce{end};
    bondN = N{end};
    bondb = b{end};

    % write to new file
    fid = fopen(bsave,'w');
    for i = 1:length(bondN)
        % Format: bondID  bondType  atom1  atom2
        fprintf(fid, '%d %d %d %f %f %d %d\n', bondAtom1(i),bondAtom2(i),bondType(i),bondLength(i),bondForce(i),bondN(i),bondb(i));
    end
    fclose(fid);

end

