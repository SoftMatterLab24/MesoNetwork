% ================= Compare2Dirs_BondDists_N_and_Lambda0.m =================
% MATLAB R2016a compatible (NO local functions in script)
%
% Reads bonds.dump (t=0) from TWO directories using your read_bonds_dump_t0.
% Plots overlay histograms (PROBABILITY per bin, not count, not density):
%   (1) Contour length normalized by b: (Lc/b) = N
%   (2) Prestretch: lambda0 = r/(N*b)
%
% Requires on path:
%   - read_bonds_dump_t0.m
%   - binned_prob_per_bin_given_edges.m

clear; clc; close all;

%% ===================== USER INPUT =====================
root = 'E:\PhD\My Research\Polydisperse_fracture\PAPER';

Samp(1) = struct( ...
    'dataDir', fullfile(root,'PD_middist_smp2_notched'), ...
    'label',   'PD\_smp2');

Samp(2) = struct( ...
    'dataDir', fullfile(root,'MD_smp2'), ...
    'label',   'Mono_smp1');

bondsDumpName = 'bonds.dump';

% Optional: override b in dump with a constant
use_b_const = false;   % set true if you want fixed b
b_const     = 1.6;

% Histogram settings
nbins_N   = 35;      % for N (=Lc/b)
nbins_lam = 35;      % for lambda0
alphaFace = 0.9;    % transparency

% Colors
col1 = [37 193 205]/255;
col2 = [255 151 151]/255;

outcol1 = [21 96 130]/255;
outcol2 = [203 54 88]/255;

% Fonts
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultFigureColor','w');

%% ===================== LOAD BOTH SAMPLES =====================
for k = 1:2
    bondsFile = fullfile(Samp(k).dataDir, bondsDumpName);
    assert(exist(bondsFile,'file')==2, 'Missing bonds.dump: %s', bondsFile);

    % Your reader: returns N_dump and b_dump if 7 columns exist
    [bond_type, i_d, j_d, r_d, f_d, N_dump, b_dump] = read_bonds_dump_t0(bondsFile); %#ok<ASGLU>

    % Require N and b in dump
    keep = isfinite(r_d) & isfinite(N_dump) & isfinite(b_dump);
    if nnz(~keep) > 0
        warning('%s: dropping %d bonds missing r/N/b in dump.', Samp(k).label, nnz(~keep));
    end

    r = r_d(keep);
    N = N_dump(keep);

    if use_b_const
        b = b_const * ones(size(N));
    else
        b = b_dump(keep);
    end

    % ---- Derived (FIXED) ----
    Lc      = N .* b;          % contour length
    Nplot   = Lc ./ b;         % = N  (this is "Lc normalized by b")
    lambda0 = r ./ (N .* b);   % = r/Lc

    Samp(k).Nplot   = Nplot(:);
    Samp(k).lambda0 = lambda0(:);

    fprintf('%s: kept %d bonds | N range [%.3g, %.3g] | lambda0 range [%.3g, %.3g]\n', ...
        Samp(k).label, numel(Nplot), min(Nplot), max(Nplot), min(lambda0), max(lambda0));
end

%% ===================== COMMON EDGES =====================
AllN   = [Samp(1).Nplot;   Samp(2).Nplot];
AllLam = [Samp(1).lambda0; Samp(2).lambda0];

AllN   = AllN(isfinite(AllN));
AllLam = AllLam(isfinite(AllLam));

tol = 1e-12;

% N edges
if (max(AllN) - min(AllN)) < tol
    N0 = AllN(1);
    edgesN = [N0-0.5, N0+0.5];    % single bin window
else
    edgesN = linspace(min(AllN), max(AllN), nbins_N+1);
end

% lambda edges
if (max(AllLam) - min(AllLam)) < tol
    L0 = AllLam(1);
    edgesL = [L0-0.05*max(L0,1), L0+0.05*max(L0,1)];
else
    edgesL = linspace(min(AllLam), max(AllLam), nbins_lam+1);
end

%% ===================== BINNED PROBABILITY (per bin) =====================
[xcN, pN1] = binned_prob_per_bin_given_edges(Samp(1).Nplot, edgesN);
[~,   pN2] = binned_prob_per_bin_given_edges(Samp(2).Nplot, edgesN);

[xcL, pL1] = binned_prob_per_bin_given_edges(Samp(1).lambda0, edgesL);
[~,   pL2] = binned_prob_per_bin_given_edges(Samp(2).lambda0, edgesL);

% ===================== ONE FIGURE, TWO SUBPLOTS: ONLY N (=Lc/b) =====================
% ===================== STACKED (2,1): ONLY bottom has x labels =====================
figure('Name','Contour length normalized by b (Lc/b = N)','Color','w');

fill1 = col1;  edge1 = outcol1;
fill2 = col2;  edge2 = outcol2;
bw = 1.00;

% optional: same y-limit for fair comparison
ymax = 1.05*max([pN1(:); pN2(:); 1e-6]);

% ---- Top: Sample 1 ----
ax1 = subplot(2,1,1); hold on; box on;

b1 = bar(xcN, pN1, bw, 'FaceColor', fill1, 'EdgeColor','none');
set(b1,'FaceAlpha',alphaFace);
o1 = bar(xcN, pN1, bw, 'FaceColor','none', 'EdgeColor',edge1, 'LineWidth',1.2);
set(o1,'HandleVisibility','off');

% ylabel('Probability per bin');
% title(sprintf('Contour length: %s', Samp(1).label));
ylim([0, ymax]);

% hide x tick labels on top subplot
set(ax1, 'XTickLabel', []);

% ---- Bottom: Sample 2 ----
ax2 = subplot(2,1,2); hold on; box on;

b2 = bar(xcN, pN2, bw, 'FaceColor', fill2, 'EdgeColor','none');
set(b2,'FaceAlpha',alphaFace);
o2 = bar(xcN, pN2, bw, 'FaceColor','none', 'EdgeColor',edge2, 'LineWidth',1.2);
set(o2,'HandleVisibility','off');

% xlabel('L_c / b \; (= N)');
% ylabel('Probability per bin');
% title(sprintf('Contour length: %s', Samp(2).label));
ylim([0, ymax]);

% ---- Make x-limits identical and keep them synced ----
xlim(ax1, [edgesN(1), edgesN(end)]);
xlim(ax2, [edgesN(1), edgesN(end)]);
linkaxes([ax1 ax2],'x');

%% ===================== PLOT 2: lambda0 OVERLAY =====================
figure('Name','Prestretch \lambda_0'); hold on; box on;

% ---- choose fill and outline colors ----
fill1 = col1;  edge1 = outcol1;     % sample 1
fill2 = col2;  edge2 = outcol2;     % sample 2

bw1 = 1.00;
bw2 = 0.85;

% filled (no edges)
b3 = bar(xcL, pL1, bw1, 'FaceColor', fill1, 'EdgeColor','none');
set(b3,'FaceAlpha',alphaFace);

b4 = bar(xcL, pL2, bw2, 'FaceColor', fill2, 'EdgeColor','none');
set(b4,'FaceAlpha',alphaFace);

% outlines (no fill)
o3 = bar(xcL, pL1, bw1, 'FaceColor','none', 'EdgeColor',edge1, 'LineWidth',1.2);
set(o3,'HandleVisibility','off');

o4 = bar(xcL, pL2, bw2, 'FaceColor','none', 'EdgeColor',edge2, 'LineWidth',1.2);
set(o4,'HandleVisibility','off');

xlabel('\lambda_0 = r/(N b)');
ylabel('Probability per bin');
title('Pre-stretch distribution (overlay)');
legend([b3 b4], Samp(1).label, Samp(2).label, 'Location','best');

if numel(edgesL)==2, xlim([edgesL(1), edgesL(2)]); end
ylim([0, 1.05*max([pL1(:); pL2(:); 1e-6])]);

%% ===================== SANITY =====================
fprintf('\nSanity: sum P(N)      sample1=%.4f, sample2=%.4f\n', sum(pN1), sum(pN2));
fprintf('Sanity: sum P(lambda) sample1=%.4f, sample2=%.4f\n', sum(pL1), sum(pL2));