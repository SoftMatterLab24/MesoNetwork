function [xc, pbin] = binned_prob_per_bin_given_edges(x, edges)
% Probability per bin: heights sum to 1, independent of bin width
[counts, ~] = histcounts(x, edges);
N_in = sum(counts);                 % # of samples that fell inside edges
if N_in == 0
    pbin = zeros(1, numel(edges)-1);
else
    pbin = counts / N_in;           % fraction in each bin
end
xc = 0.5*(edges(1:end-1) + edges(2:end));
end
