function [id, iatom, jatom, N, b] = parse_bond_table(tfile)

%% Parsing bond table

fid = fopen(tfile);

% Skip first five lines
for skip = 1:5
    tline = fgetl(fid);
end
 
% Read bond data
C = textscan(fid,'%u %u %u %f %f');

id      = C{1};
iatom   = C{2};
jatom   = C{3};
N       = C{4};
b       = C{5};


fclose(fid);

end