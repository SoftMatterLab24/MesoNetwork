% PlotBondDistributionsMulti_Lc.m
% MATLAB R2016a compatible
% Multi-sample equilibrium statistics with consistent binning & averaging
%
% Reads each sample:
%   bonds.dump  (columns: bond_type, i, j, r, f [, b1=N_kuhn, b2=b_kuhn])
%   bond.table  (columns: index, i, j, N_kuhn, b_kuhn)  -- only if needed
%
% Outputs (all bonds after matching N,b when needed):
%   - End-to-end length r distribution                (continuous)
%   - Contour length Lc = N * b distribution          (continuous)   <-- NEW
%   - Pre-stretch lambda0 = r/(N*b) distribution      (continuous)
%
% Averaging modes:
%   'pooled'           : pool all bonds across samples, then histogram
%   'per-sample-mean'  : histogram each sample on common bins, then average PDFs

clear; clc; close all;

%% ===================== USER KNOBS =====================
% --- Samples list (add as many as you want) ---
Samples = struct( ...
  'dataDir'      , {}, ...
  'bondsDumpFile', {}, ...
  'bondTableFile', {}, ...
  'has_b_in_dump', {});     % true: bonds.dump has b1,b2; false: read bond.table

% EXAMPLE: add samples (edit paths as needed)
root = 'E:\PhD\My Research\Polydisperse_fracture\PAPER';
% Samples(end+1) = struct( ...
%   'dataDir',       fullfile(root,'bimodal_equilibrated_samples'), ...
%   'bondsDumpFile', fullfile(root,'bimodal_equilibrated_samples','bimodal_bondequib_001.dump'), ...
%   'bondTableFile', fullfile(root,'bimodal_equilibrated_samples','bond.table'), ...
%   'has_b_in_dump', true);
% 
Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'PD_smp2'), ...
  'bondsDumpFile', fullfile(root,'PD_smp2','bonds.dump'), ...
  'bondTableFile', fullfile(root,'PD_smp2','bond.table'), ...
  'has_b_in_dump', true);
% % 
% Samples(end+1) = struct( ...
%   'dataDir',       fullfile(root,'MD_smp4'), ...
%   'bondsDumpFile', fullfile(root,'MD_smp4','bonds.dump'), ...
%   'bondTableFile', fullfile(root,'MD_smp4','bond.table'), ...
%   'has_b_in_dump', true);
% % 
% Samples(end+1) = struct( ...
%   'dataDir',       fullfile(root,'MD_smp5'), ...
%   'bondsDumpFile', fullfile(root,'MD_smp5','bonds.dump'), ...
%   'bondTableFile', fullfile(root,'MD_smp5','bond.table'), ...
%   'has_b_in_dump', true);

% --- Plot look ---
plotMode        = 'hist';            % 'hist' or 'line'
avg_mode        = 'pooled';          % 'pooled' or 'per-sample-mean'

% --- Consistent ranges (leave [] to auto from all samples) ---
r_range         = [];                % e.g., [0, 250]
Lc_range        = [];                % e.g., [0, 200]   (NEW)
lambda_range    = [];                % e.g., [0, 2]

% --- Unified bin count for ALL THREE plots ---
nbins_common    = 35;                % one knob to rule them all

% --- Kuhn length handling ---
use_b_const     = true;              % true: override/ignore per-bond b and use b_const
b_const         = 1.6;               % constant Kuhn length (nm, or your units)

% --- Misc ---
show_per_sample_overlays = false;
col_main  = [221 132 146]/255;
color_line = [203 54 88]/255;

%% ===================== LOAD ALL SAMPLES (PASS 1) =====================
ns = numel(Samples);
assert(ns>=1, 'Please add at least one Sample in the "Samples" struct array.');

All_r = []; All_N = []; All_b = []; All_Lc = []; All_lambda = [];
PerSample = struct('r',[],'N',[],'b',[],'Lc',[],'lambda0',[],'bond_type',[]);
All_bond_type = [];

for s = 1:ns
    fprintf('\n=== Sample %d/%d ===\n', s, ns);
    S = Samples(s);
    assert(exist(S.bondsDumpFile,'file')==2, 'Missing bonds.dump: %s', S.bondsDumpFile);
    if ~S.has_b_in_dump
        assert(exist(S.bondTableFile,'file')==2, 'Missing bond.table (needed): %s', S.bondTableFile);
    end

    % ---- Read t=0 dump block ----
    [bond_type, i_d, j_d, r_d, f_d, N_in_dump, b_in_dump] = read_bonds_dump_t0(S.bondsDumpFile);

    % ---- Get N, b per bond ----
    if S.has_b_in_dump
        keep = ~isnan(N_in_dump) & ~isnan(b_in_dump);
        if ~all(keep)
            warning('Sample %d: %d / %d bonds missing N,b in dump; dropping those.', s, nnz(~keep), numel(keep));
        end
        r   = r_d(keep);
        N   = N_in_dump(keep);
        if use_b_const
            b = b_const * ones(size(N));
        else
            b = b_in_dump(keep);
        end
        btp = bond_type(keep);
    else
        [N_tab, b_tab] = read_bond_table_map(S.bondTableFile);
        [r, N, b, btp] = match_dump_with_table(i_d, j_d, r_d, bond_type, N_tab, b_tab);
        if use_b_const
            b = b_const * ones(size(b));
        end
    end

    % ---- Derived quantities ----
    Lc = N .* b;                 % contour length (continuous)
    lambda0 = r ./ Lc;           % pre-stretch

    % ---- Store per-sample ----
    PerSample(s).r        = r;
    PerSample(s).N        = N;
    PerSample(s).b        = b;
    PerSample(s).Lc       = Lc;
    PerSample(s).lambda0  = lambda0;
    PerSample(s).bond_type= btp;

    % ---- Accumulate globals for auto-ranges ----
    All_r      = [All_r; r];        %#ok<AGROW>
    All_N      = [All_N; N];        %#ok<AGROW>
    All_b      = [All_b; b];        %#ok<AGROW>
    All_Lc     = [All_Lc; Lc];      %#ok<AGROW>
    All_lambda = [All_lambda; lambda0]; %#ok<AGROW>
    
    
    All_bond_type = [All_bond_type; btp]; %#ok<AGROW>
    fprintf('Sample %d: kept %d bonds. r[min,mean,max]=[%.3g, %.3g, %.3g], Lc[min,max]=[%.3g,%-.3g], lambda[min,max]=[%.3g, %.3g]\n', ...
        s, numel(r), min(r), mean(r), max(r), min(Lc), max(Lc), min(lambda0), max(lambda0));
end

%% ===================== BUILD COMMON BINS =====================
% r (continuous)
if isempty(r_range)
    r_range = [min(All_r), max(All_r)];
end
edges_r = linspace(r_range(1), r_range(2), nbins_common+1);
xc_r    = 0.5*(edges_r(1:end-1) + edges_r(2:end));   % centers

% Lc (continuous)  <-- NEW
if isempty(Lc_range)
    Lc_range = [min(All_Lc), max(All_Lc)];
end
edges_Lc = linspace(Lc_range(1), Lc_range(2), nbins_common+1);
xc_Lc    = 0.5*(edges_Lc(1:end-1) + edges_Lc(2:end));

% lambda (continuous)
if isempty(lambda_range)
    lambda_range = [min(All_lambda), max(All_lambda)];
end
edges_l = linspace(lambda_range(1), lambda_range(2), nbins_common+1);
xc_l    = 0.5*(edges_l(1:end-1) + edges_l(2:end));

%% ===================== HISTOGRAMS / PDFs =====================
% Precompute pooled vectors (used for "pooled" avg and for type-split)
r_pool   = concatenate_field(PerSample,'r');
Lc_pool  = concatenate_field(PerSample,'Lc');
l_pool   = concatenate_field(PerSample,'lambda0');
bt_pool  = concatenate_field(PerSample,'bond_type');
unique_types = unique(All_bond_type);

switch avg_mode
    case 'pooled'
        [~, prob_r ] = binned_prob_per_bin_given_edges(r_pool , edges_r );
        [~, prob_Lc] = binned_prob_per_bin_given_edges(Lc_pool, edges_Lc);
        [~, prob_l ] = binned_prob_per_bin_given_edges(l_pool , edges_l );

        y_r_mean  = prob_r;
        y_Lc_mean = prob_Lc;
        y_l_mean  = prob_l;

    case 'per-sample-mean'
        Y_r  = zeros(ns, numel(xc_r ));
        Y_Lc = zeros(ns, numel(xc_Lc));
        Y_l  = zeros(ns, numel(xc_l ));
        for s = 1:ns
            [~, pr ] = binned_prob_per_bin_given_edges(PerSample(s).r      , edges_r );
            [~, pLc] = binned_prob_per_bin_given_edges(PerSample(s).Lc     , edges_Lc);
            [~, pl ] = binned_prob_per_bin_given_edges(PerSample(s).lambda0, edges_l );
            Y_r(s,:)  = pr;
            Y_Lc(s,:) = pLc;
            Y_l(s,:)  = pl;
        end
        y_r_mean  = mean(Y_r ,1);
        y_Lc_mean = mean(Y_Lc,1);
        y_l_mean  = mean(Y_l ,1);

    otherwise
        error('Unknown avg_mode: %s', avg_mode);
end

% ---------- OPTIONAL: split by bond type (up to 2 types) ----------
y_r_t1 = []; y_r_t2 = [];
y_Lc_t1 = []; y_Lc_t2 = [];
y_l_t1 = []; y_l_t2 = [];

if numel(unique_types) == 2
    t1 = unique_types(1);
    t2 = unique_types(2);

    % Bin widths (constant, but keep general)
    dr  = diff(edges_r);
    dLc = diff(edges_Lc);
    dl  = diff(edges_l);

    switch avg_mode
        case 'pooled'
            % ------- POOL ALL BONDS, NORMALIZE BY TOTAL CHAIN COUNT -------
            totalBonds = numel(r_pool);

            idx1 = (bt_pool == t1);
            idx2 = (bt_pool == t2);

            r1  = r_pool(idx1);  r2  = r_pool(idx2);
            Lc1 = Lc_pool(idx1); Lc2 = Lc_pool(idx2);
            l1  = l_pool(idx1);  l2  = l_pool(idx2);

            % r
            counts_r1 = histcounts(r1 , edges_r);
            counts_r2 = histcounts(r2 , edges_r);
            y_r_t1 = counts_r1 / totalBonds ./ dr;
            y_r_t2 = counts_r2 / totalBonds ./ dr;

            % Lc
            counts_Lc1 = histcounts(Lc1, edges_Lc);
            counts_Lc2 = histcounts(Lc2, edges_Lc);
            y_Lc_t1 = counts_Lc1 / totalBonds ./ dLc;
            y_Lc_t2 = counts_Lc2 / totalBonds ./ dLc;

            % lambda0
            counts_l1 = histcounts(l1 , edges_l);
            counts_l2 = histcounts(l2 , edges_l);
            y_l_t1 = counts_l1 / totalBonds ./ dl;
            y_l_t2 = counts_l2 / totalBonds ./ dl;

        case 'per-sample-mean'
            % ------- PER-SAMPLE PDF, STILL NORMALIZED BY TOTAL CHAINS IN SAMPLE -------
            nb_r  = numel(xc_r );
            nb_Lc = numel(xc_Lc);
            nb_l  = numel(xc_l );

            Y_r_t1  = zeros(ns, nb_r );
            Y_r_t2  = zeros(ns, nb_r );
            Y_Lc_t1 = zeros(ns, nb_Lc);
            Y_Lc_t2 = zeros(ns, nb_Lc);
            Y_l_t1  = zeros(ns, nb_l );
            Y_l_t2  = zeros(ns, nb_l );

            for s = 1:ns
                bt_s = PerSample(s).bond_type;
                rs   = PerSample(s).r;
                Lcs  = PerSample(s).Lc;
                ls   = PerSample(s).lambda0;

                Ns = numel(rs); % total chains in this sample

                idx1 = (bt_s == t1);
                idx2 = (bt_s == t2);

                % --- type 1 ---
                if any(idx1)
                    c_r1  = histcounts(rs(idx1) , edges_r );
                    c_Lc1 = histcounts(Lcs(idx1), edges_Lc);
                    c_l1  = histcounts(ls(idx1) , edges_l );
                    pr1   = c_r1  / Ns ./ dr;
                    pLc1  = c_Lc1 / Ns ./ dLc;
                    pl1   = c_l1  / Ns ./ dl;
                else
                    pr1  = zeros(1, nb_r );
                    pLc1 = zeros(1, nb_Lc);
                    pl1  = zeros(1, nb_l );
                end

                % --- type 2 ---
                if any(idx2)
                    c_r2  = histcounts(rs(idx2) , edges_r );
                    c_Lc2 = histcounts(Lcs(idx2), edges_Lc);
                    c_l2  = histcounts(ls(idx2) , edges_l );
                    pr2   = c_r2  / Ns ./ dr;
                    pLc2  = c_Lc2 / Ns ./ dLc;
                    pl2   = c_l2  / Ns ./ dl;
                else
                    pr2  = zeros(1, nb_r );
                    pLc2 = zeros(1, nb_Lc);
                    pl2  = zeros(1, nb_l );
                end

                Y_r_t1(s,:)  = pr1;
                Y_r_t2(s,:)  = pr2;
                Y_Lc_t1(s,:) = pLc1;
                Y_Lc_t2(s,:) = pLc2;
                Y_l_t1(s,:)  = pl1;
                Y_l_t2(s,:)  = pl2;
            end

            % Mean across samples
            y_r_t1  = mean(Y_r_t1 ,1);
            y_r_t2  = mean(Y_r_t2 ,1);
            y_Lc_t1 = mean(Y_Lc_t1,1);
            y_Lc_t2 = mean(Y_Lc_t2,1);
            y_l_t1  = mean(Y_l_t1 ,1);
            y_l_t2  = mean(Y_l_t2 ,1);
    end
end

%% ===================== PLOTS =====================
figsize = [1 1 3 3];
set(0,'DefaultFigureUnits','inches','DefaultFigurePosition',figsize);
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',15);

% ---- r ----
figure('Name','End-to-end length r (averaged)'); hold on;
set(gca,'FontName','Times New Roman','FontSize',15);

if numel(unique_types) == 2
    t1 = unique_types(1);
    t2 = unique_types(2);
    col1 = [63 61 153]/255;
    col2 = [100 153 202]/255;

    if strcmp(plotMode,'hist')
        bar(xc_r, y_r_t1, 'FaceColor',col1, 'EdgeColor','none', 'BarWidth',1);
        bar(xc_r, y_r_t2, 'FaceColor',col2, 'EdgeColor','none', 'BarWidth',1, 'FaceAlpha',0.6);
    else
        plot(xc_r, y_r_t1, 'LineWidth',1.8, 'Color',col1);
        plot(xc_r, y_r_t2, 'LineWidth',1.8, 'Color',col2);
    end
    legend(sprintf('type %d',t1), sprintf('type %d',t2), 'Location','Best');
else
    if strcmp(plotMode,'hist')
        bar(xc_r, y_r_mean, 'FaceColor',col_main, 'EdgeColor',color_line, 'BarWidth',1);
    else
        plot(xc_r, y_r_mean, 'LineWidth',1.8, 'Color', color_line);
    end
end

xlabel('End-to-end length r'); ylabel('Probability per bin');
title(sprintf('r distribution (%s)', avg_mode)); box on; grid off;
axis([10 100 0 0.4])

% ---- Lc = N*b (monodisperse-aware) ----

all_N_unique = unique(All_N);

figure('Name','Contour length L_c = N \times b (averaged)'); hold on;
set(gca,'FontName','Times New Roman','FontSize',15);

% Pure monodisperse (single N) ? keep the "clean single bar" behavior
if numel(all_N_unique) == 1 && numel(unique_types) == 1
    Lc0 = all_N_unique * b_const;  % since you're using b_const

    bar(Lc0, 1, 'FaceColor', col_main, 'EdgeColor', color_line, 'BarWidth', 1);

    xlabel('Contour length L_c = N \times b');
    ylabel('Probability');
    title('L_c distribution (monodisperse)');
    xlim([Lc0 - 0.5*Lc0, Lc0 + 0.5*Lc0]);
    ylim([0, 1.1]);
    box on; grid off;

    fprintf('Detected monodisperse N = %.4g -> L_c = %.4g. Plotted single bar.\n', ...
        all_N_unique, Lc0);

else
    if numel(unique_types) == 2
        t1 = unique_types(1);
        t2 = unique_types(2);
        col1 = [63 61 153]/255;
    col2 = [100 153 202]/255;

        if strcmp(plotMode,'hist')
            bar(xc_Lc, y_Lc_t1, 'FaceColor',col1, 'EdgeColor','none', 'BarWidth',1);
            bar(xc_Lc, y_Lc_t2, 'FaceColor',col2, 'EdgeColor','none', 'BarWidth',1, 'FaceAlpha',0.6);
        else
            plot(xc_Lc, y_Lc_t1, 'LineWidth',1.8, 'Color',col1);
            plot(xc_Lc, y_Lc_t2, 'LineWidth',1.8, 'Color',col2);
        end
        legend(sprintf('type %d',t1), sprintf('type %d',t2), 'Location','Best');
    else
        if strcmp(plotMode,'hist')
            bar(xc_Lc, y_Lc_mean, 'FaceColor',col_main, 'EdgeColor',color_line, 'BarWidth',1);
        else
            plot(xc_Lc, y_Lc_mean, 'LineWidth',1.8, 'Color', color_line);
        end
    end

    xlabel('Contour length L_c = N \times b');
    ylabel('Probability per bin');
    title(sprintf('L_c distribution (%s)', avg_mode));
    box on; grid off;
end


% ---- lambda0 ----
figure('Name','Pre-stretch \lambda_0 (averaged)'); hold on;
set(gca,'FontName','Times New Roman','FontSize',15);

if numel(unique_types) == 2
    t1 = unique_types(1);
    t2 = unique_types(2);
    col1 = [63 61 153]/255;
    col2 = [100 153 202]/255;

    if strcmp(plotMode,'hist')
        bar(xc_l, y_l_t1, 'FaceColor',col1, 'EdgeColor','none', 'BarWidth',1);
        bar(xc_l, y_l_t2, 'FaceColor',col2, 'EdgeColor','none', 'BarWidth',1, 'FaceAlpha',0.6);
    else
        plot(xc_l, y_l_t1, 'LineWidth',1.8, 'Color',col1);
        plot(xc_l, y_l_t2, 'LineWidth',1.8, 'Color',col2);
    end
    legend(sprintf('type %d',t1), sprintf('type %d',t2), 'Location','Best');
else
    if strcmp(plotMode,'hist')
        bar(xc_l, y_l_mean, 'FaceColor',col_main, 'EdgeColor',color_line, 'BarWidth',1);
    else
        plot(xc_l, y_l_mean, 'LineWidth',1.8, 'Color', color_line);
    end
end

xlabel('\lambda_0 = r / (N \cdot b)'); ylabel('Probability per bin');
title(sprintf('Pre-stretch distribution (%s)', avg_mode)); box on; grid off;

%% ===================== 2D HISTOGRAMS (pooled) =====================

% --- custom grey ? red colormap ---
ncol = 256;
c_lo = [0.9 0.9 0.9];   % light grey
c_hi = [1.0 0.0 0.0];   % red
cmap_gr2red = [linspace(c_lo(1),c_hi(1),ncol)', ...
               linspace(c_lo(2),c_hi(2),ncol)', ...
               linspace(c_lo(3),c_hi(3),ncol)'];

% ---------- 2D hist: lambda (x) vs r (y) ----------
figure('Name','2D hist: \lambda_0 vs r','Color','w');
h1 = histogram2(All_lambda, All_r, edges_l, edges_r);  % uses your common bins
% h1.DisplayStyle = 'tile';
% h1.ShowEmptyBins = 'on';
h1.Normalization = 'probability';   % probability intensity
set(gca,'FontName','Times New Roman','FontSize',15);
xlabel('\lambda_0 = r / (N \cdot b)');
ylabel('End-to-end length r');
title('P(r,\lambda_0)');
colormap(gca, cmap_gr2red);
colorbar;

% ---------- 2D hist: lambda (x) vs L_c (y) ----------
figure('Name','2D hist: \lambda_0 vs L_c','Color','w');
h2 = histogram2(All_lambda, All_Lc, edges_l, edges_Lc);
% h2.DisplayStyle = 'tile';
% h2.ShowEmptyBins = 'on';
h2.Normalization = 'probability';
set(gca,'FontName','Times New Roman','FontSize',15);
xlabel('\lambda_0 = r / (N \cdot b)');
ylabel('Contour length L_c = N \cdot b');
title('P(L_c,\lambda_0)');
colormap(gca, cmap_gr2red);
colorbar;


%% ===================== TEXT SUMMARY =====================
fprintf('\n=== GLOBAL SUMMARY (using common ranges) ===\n');
fprintf('r   range used: [%.3g, %.3g]   (nbins=%d)\n', r_range(1), r_range(2), nbins_common);
fprintf('L_c range used: [%.3g, %.3g]   (nbins=%d)\n', Lc_range(1), Lc_range(2), nbins_common);
fprintf('lam range used: [%.3g, %.3g]   (nbins=%d)\n', lambda_range(1), lambda_range(2), nbins_common);
fprintf('Averaging mode: %s | plotMode: %s | b_const=%g (use_b_const=%d)\n', avg_mode, plotMode, b_const, use_b_const);
fprintf('sum(prob_r)=%.6f, sum(prob_Lc)=%.6f, sum(prob_l)=%.6f\n', sum(y_r_mean), sum(y_Lc_mean), sum(y_l_mean));
