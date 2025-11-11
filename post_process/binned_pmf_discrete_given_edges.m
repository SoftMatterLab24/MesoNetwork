function [xc, pmf] = binned_pmf_discrete_given_edges(x, edges)
% Integer-centered PMF (probabilities that sum to 1)
counts = histcounts(x, edges);
pmf = counts / max(1,sum(counts));
xc  = edges(1:end-1)+0.5; % integer centers
end