function [r, N, b, btp] = match_dump_with_table(i_d, j_d, r_d, bond_type, MapN, Mapb)
n = numel(r_d);
N = nan(n,1); b = nan(n,1);
for k=1:n
    key = sprintf('%d_%d', min(i_d(k),j_d(k)), max(i_d(k),j_d(k)));
    if MapN.isKey(key) && Mapb.isKey(key)
        N(k) = MapN(key); b(k)=Mapb(key);
    end
end
keep = ~isnan(N) & ~isnan(b);
if ~all(keep)
    warning('Dropped %d / %d bonds without N,b match.', nnz(~keep), n);
end
r   = r_d(keep);
N   = N(keep);
b   = b(keep);
btp = bond_type(keep);
end