function [xc, pmf_vals] = binned_pmf_discrete(x, edges)
% For integer-like data (e.g., N): PMF line at integer centers from given edges.
x = x(~isnan(x) & isfinite(x));
cnts = histcounts(x, edges);
pmf_vals = cnts / sum(cnts);
xc = 0.5*(edges(1:end-1) + edges(2:end)); % integer centers (Nmin:Nmax)
end