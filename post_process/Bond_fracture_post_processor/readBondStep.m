function [ok, timestep, Nbonds, xlo, xhi, ylo, yhi, zlo, zhi, ...
          typeB, id1B, id2B, lenB, forceB, NB, bB] = readBondStep(fid)

    ok = false;
    timestep = []; Nbonds = [];
    xlo = []; xhi = []; ylo = []; yhi = []; zlo = []; zhi = [];
    typeB = []; id1B = []; id2B = [];
    lenB = []; forceB = []; NB = []; bB = [];

    % ----- find "ITEM: TIMESTEP" -----
    line = fgetl(fid);
    if ~ischar(line)
        return; % EOF
    end
    % If this line isn't the TIMESTEP header, search forward
    while ischar(line) && ~strncmp(line, 'ITEM: TIMESTEP', 14)
        line = fgetl(fid);
        if ~ischar(line)
            return; % EOF
        end
    end

    % ----- timestep value -----
    line = fgetl(fid);
    if ~ischar(line), return; end
    timestep = sscanf(line, '%d', 1);

    % ----- "ITEM: NUMBER OF ENTRIES" -----
    line = fgetl(fid);
    if ~ischar(line), return; end
    % (we don't care what exactly it says, just skip it)

    % ----- number of bonds -----
    line = fgetl(fid);
    if ~ischar(line), return; end
    Nbonds = sscanf(line, '%d', 1);

    % ----- find "ITEM: BOX BOUNDS" -----
    line = fgetl(fid);
    if ~ischar(line), return; end
    while ischar(line) && ~strncmp(line, 'ITEM: BOX BOUNDS', 16)
        line = fgetl(fid);
        if ~ischar(line)
            return; % EOF
        end
    end

    % ----- read x/y/z bounds -----
    [xlo, xhi] = readBoundsLine(fid);
    [ylo, yhi] = readBoundsLine(fid);
    [zlo, zhi] = readBoundsLine(fid);

    % ----- "ITEM: ENTRIES ..." -----
    line = fgetl(fid);
    if ~ischar(line), return; end
    % expected: ITEM: ENTRIES c_1[1] c_1[2] ...

    % ----- Nbonds lines of numeric data -----
    % type id1 id2 length force N_kuhn b_kuhn
    fmt  = '%f %f %f %f %f %f %f';
    data = textscan(fid, fmt, Nbonds, ...
        'Delimiter', ' ', 'CollectOutput', true);

    data = data{1};
    if size(data,1) ~= Nbonds
        error('Failed to read %d bond lines at timestep %d', Nbonds, timestep);
    end

    typeB  = data(:,1);
    id1B   = data(:,2);
    id2B   = data(:,3);
    lenB   = data(:,4);
    forceB = data(:,5);
    NB     = data(:,6);
    bB     = data(:,7);

    ok = true;
end
