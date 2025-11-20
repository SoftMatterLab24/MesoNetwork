function [x,y,id,ty,mol,natoms,bAtom1,bAtom2,bType,bLength,bForce,xlims,ylims,zlims] = parse_equidump_full_fun(afile,bfile)

%afile = "atomDump.dump";
%bfile = "bondsDump.dump";

%% Parsing atoms

fid = fopen(afile);

ii = 0;
while 1 == 1

    ii = ii + 1;

    % Skip first line
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end

    % Timestep
    tstep(ii) = str2double(fgetl(fid));

    % Skip
    tline = fgetl(fid);

    % Number of atoms in the system
    natoms(ii) = str2double(fgetl(fid));

    % Skip header
    tline = fgetl(fid);

    % Simulation dimensions
    xlims{ii} = str2double(strsplit(fgetl(fid)));
    ylims{ii} = str2double(strsplit(fgetl(fid)));
    zlims{ii} = str2double(strsplit(fgetl(fid)));

    % Column legend
    columns = fgetl(fid);

    % Read the rest of the file 
    C = textscan(fid,'%u %u %u %f %f %f %f %f %f');

    id{ii}  = C{1};
    mol{ii} = C{2};
    ty{ii}  = C{3};
    x{ii}   = C{4};
    y{ii}   = C{5};
    %z{ii}   = 0;

end

fclose(fid);

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
    b1{ii} = C{6};
    b2{ii} = C{7};
    % bForcez{ii} = C{8};
    % bEngpot{ii} = C{6};


end
fclose(fid);


