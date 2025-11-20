function [ok, timestep, Natoms, xlo, xhi, ylo, yhi, zlo, zhi, ...
          atomID, xs, ys, zs, fx, fy, fz] = readAtomStep(fid)

    ok = false;
    timestep = []; Natoms = [];
    xlo = []; xhi = []; ylo = []; yhi = []; zlo = []; zhi = [];
    atomID = []; xs = []; ys = []; zs = [];
    fx = []; fy = []; fz = [];

    % ----- find "ITEM: TIMESTEP" -----
    line = fgetl(fid);
    if ~ischar(line)
        return; % EOF
    end
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

    % ----- "ITEM: NUMBER OF ATOMS" -----
    line = fgetl(fid);
    if ~ischar(line), return; end

    % ----- Natoms -----
    line = fgetl(fid);
    if ~ischar(line), return; end
    Natoms = sscanf(line, '%d', 1);

    % ----- find "ITEM: BOX BOUNDS" -----
    line = fgetl(fid);
    if ~ischar(line), return; end
    while ischar(line) && ~strncmp(line, 'ITEM: BOX BOUNDS', 16)
        line = fgetl(fid);
        if ~ischar(line)
            return; % EOF
        end
    end

    % ----- x/y/z bounds -----
    [xlo, xhi] = readBoundsLine(fid);
    [ylo, yhi] = readBoundsLine(fid);
    [zlo, zhi] = readBoundsLine(fid);

    % ----- "ITEM: ATOMS ..." -----
    line = fgetl(fid);
    if ~ischar(line), return; end
    % expected: ITEM: ATOMS id mol type xs ys zs fx fy fz

    % ----- Natoms lines of numeric data -----
    % id mol type xs ys zs fx fy fz
    fmt  = '%f %f %f %f %f %f %f %f %f';
    data = textscan(fid, fmt, Natoms, ...
        'Delimiter', ' ', 'CollectOutput', true);

    data = data{1};
    if size(data,1) ~= Natoms
        error('Failed to read %d atom lines at timestep %d', Natoms, timestep);
    end

    atomID = data(:,1);
    xs     = data(:,4);
    ys     = data(:,5);
    zs     = data(:,6);
    fx     = data(:,7);
    fy     = data(:,8);
    fz     = data(:,9);

    ok = true;
end
