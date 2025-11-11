function [xc, pdf] = binned_pdf_cont_given_edges(x, edges)
% Returns a probability density (? pdf dx = 1)

counts = histcounts(x, edges);
dx = diff(edges);
N = sum(counts);
pdf = counts ./ (N .* dx);
xc  = 0.5*(edges(1:end-1) + edges(2:end));

% guard against NaN or Inf
pdf(~isfinite(pdf)) = 0;
end
