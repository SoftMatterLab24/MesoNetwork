function [MapN, Mapb] = read_bond_table_map(fname)
% Returns containers.Map from 'imin_imax' -> N or b
fid = fopen(fname,'r'); assert(fid>0,'Could not open %s',fname);
i_t=[]; j_t=[]; N_t=[]; b_t=[];
ln = fgetl(fid);
while ischar(ln)
    s = strtrim(ln);
    if isempty(s) || s(1)=='#' || strncmpi(s,'KEY',3) || s(1)=='N'
        ln = fgetl(fid); continue;
    end
    vals = sscanf(s,'%f %f %f %f %f');
    if numel(vals)==5
        i_t(end+1,1)=vals(2); j_t(end+1,1)=vals(3);
        N_t(end+1,1)=vals(4); b_t(end+1,1)=vals(5);
    end
    ln = fgetl(fid);
end
fclose(fid);
MapN = containers.Map('KeyType','char','ValueType','double');
Mapb = containers.Map('KeyType','char','ValueType','double');
for k=1:numel(i_t)
    key = sprintf('%d_%d', min(i_t(k),j_t(k)), max(i_t(k),j_t(k)));
    MapN(key) = N_t(k);
    Mapb(key) = b_t(k);
end
end