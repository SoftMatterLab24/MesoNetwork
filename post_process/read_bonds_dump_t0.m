function [bond_type, i_d, j_d, r_d, f_d, N_dump, b_dump] = read_bonds_dump_t0(fname)
% Parse only t=0 block from a LAMMPS-style bond dump, OR a plain numeric file.
% Accepts 5 or 7 numeric columns per entry:
%   5 cols: [type, i, j, r, f]
%   7 cols: [type, i, j, r, f, b1(=N), b2(=b)]

fid = fopen(fname,'r');
assert(fid>0,'Cannot open %s',fname);

% -------- Detect if file has LAMMPS "ITEM:" headers or is plain numeric -----
firstLine = '';
while ~feof(fid)
    firstLine = fgetl(fid);
    if ~ischar(firstLine), break; end
    if ~isempty(strtrim(firstLine))
        break; % first non-empty line
    end
end

hasHeader = false;
if ischar(firstLine)
    if strncmp(strtrim(firstLine),'ITEM: TIMESTEP',14)
        hasHeader = true;
    end
end

% Rewind to start for actual parsing
fseek(fid,0,'bof');

B      = [];
N_dump = [];
b_dump = [];

if hasHeader
    % ================= LAMMPS DUMP MODE (original logic) ==================
    reading = false;
    done    = false;
    while ~feof(fid) && ~done
        ln = fgetl(fid);
        if ~ischar(ln), break; end

        if strncmp(ln,'ITEM: TIMESTEP',14)
            ts = str2double(fgetl(fid));
            reading = (ts==0);
            continue;
        end

        if reading && strncmp(ln,'ITEM: ENTRIES',13)
            while true
                pos = ftell(fid);
                L = fgetl(fid);
                if ~ischar(L) || strncmp(L,'ITEM:',5)
                    fseek(fid,pos,'bof');
                    break;
                end
                vals = sscanf(L,'%f');
                if numel(vals)==5 || numel(vals)==7
                    B(end+1,1:5) = vals(1:5).'; %#ok<SAGROW>
                    if numel(vals)==7
                        N_dump(end+1,1) = vals(6);
                        b_dump(end+1,1) = vals(7);
                    else
                        N_dump(end+1,1) = NaN;
                        b_dump(end+1,1) = NaN;
                    end
                end
            end
            done = true;
        end
    end

else
    % ================= PLAIN NUMERIC FILE MODE ==================
    % No headers, just rows of numbers. Treat everything as t = 0.
    while ~feof(fid)
        L = fgetl(fid);
        if ~ischar(L), break; end
        if isempty(strtrim(L)), continue; end

        vals = sscanf(L,'%f');
        if numel(vals)==5 || numel(vals)==7
            B(end+1,1:5) = vals(1:5).'; %#ok<SAGROW>
            if numel(vals)==7
                N_dump(end+1,1) = vals(6);
                b_dump(end+1,1) = vals(7);
            else
                N_dump(end+1,1) = NaN;
                b_dump(end+1,1) = NaN;
            end
        end
    end
end

fclose(fid);
assert(~isempty(B),'No numeric entries read from timestep 0 (or file) in %s', fname);

bond_type = B(:,3);
i_d       = B(:,1);
j_d       = B(:,2);
r_d       = B(:,4);
f_d       = B(:,5);
end
