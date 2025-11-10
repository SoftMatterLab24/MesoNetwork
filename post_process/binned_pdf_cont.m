function [xc, pdf_vals] = binned_pdf_cont(x, nbins)
% For continuous data: make a PDF line from counts using uniform bins.
% Returns bin centers xc and density pdf_vals so that sum(pdf*binwidth)?1.
x = x(~isnan(x) & isfinite(x));
xmin = min(x); xmax = max(x);
if xmin==xmax
    xc = xmin; pdf_vals = 1; return;
end
edges = linspace(xmin, xmax, nbins+1);
binw  = edges(2) - edges(1);
cnts  = histcounts(x, edges);               % counts per bin
pdf_vals = cnts / (numel(x) * binw);        % convert to PDF
xc = 0.5*(edges(1:end-1) + edges(end:-1:2)); % WRONG; fix to:
xc = 0.5*(edges(1:end-1) + edges(2:end));   % bin centers
end