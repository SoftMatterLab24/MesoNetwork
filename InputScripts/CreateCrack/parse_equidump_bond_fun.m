function [bAtom1,bAtom2,bType,bLength,bForce,N,b] = parse_equidump_bond_fun(bfile)

%afile = "atomDump.dump";
%bfile = "bondsDump.dump";

%% Parsing bonds

fid = fopen(bfile);

ii = 0;
while 1 == 1

    ii = ii + 1;

    % Skip first line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end

    % Timestep
    tbstep(ii) = str2double(fgetl(fid));

    % Skip
    tline = fgetl(fid);

    % Number of atoms in the system
    nbonds(ii) = str2double(fgetl(fid));
    if nbonds(ii) == 0
        % Skip
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        tline = fgetl(fid);
        continue
    end

    % Skip header
    tline = fgetl(fid);

    % Simulation dimensions
    xlimsb{ii} = str2double(strsplit(fgetl(fid)));
    ylimsb{ii} = str2double(strsplit(fgetl(fid)));
    zlimsb{ii} = str2double(strsplit(fgetl(fid)));

    % Column legend
    columns = fgetl(fid);

    % Read the rest of the file 
    C = textscan(fid,'%u %u %u %f %f %d %d');

    bType{ii}   = C{1};
    bAtom1{ii}  = C{2};
    bAtom2{ii}  = C{3};
    bLength{ii} = C{4};
    bForce{ii}  = C{5};
    N{ii}       = C{6};
    b{ii}       = C{7};
    %b1{ii} = C{6};
    %b2{ii} = C{7};
    % bForcez{ii} = C{8};
    % bEngpot{ii} = C{6};


end
fclose(fid);


