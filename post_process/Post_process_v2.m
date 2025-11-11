% PlotBondDistributionsMulti.m
% MATLAB R2016a compatible
% Multi-sample equilibrium statistics with consistent binning & averaging
%
% Reads each sample:
%   bonds.dump  (columns: bond_type, i, j, r, f [, b1=N_kuhn, b2=b_kuhn])
%   bond.table  (columns: index, i, j, N_kuhn, b_kuhn)  -- only if needed
%
% Outputs (all bonds after matching N,b when needed):
%   - End-to-end length r distribution
%   - Contour length (Kuhn segments) N distribution
%   - Pre-stretch lambda0 = r/(N*b) distribution
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
Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'sample_alpine_huge_erate_1e-6'), ...
  'bondsDumpFile', fullfile(root,'sample_alpine_huge_erate_1e-6','bonds.dump'), ...
  'bondTableFile', fullfile(root,'sample_alpine_huge_erate_1e-6','bond.table'), ...
  'has_b_in_dump', false);   % this one needs bond.table for N,b

Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'PD_smp2'), ...
  'bondsDumpFile', fullfile(root,'PD_smp2','bonds.dump'), ...
  'bondTableFile', fullfile(root,'PD_smp2','bond.table'), ...
  'has_b_in_dump', true);   % this one does not need bond.table for N,b

Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'PD_smp3'), ...
  'bondsDumpFile', fullfile(root,'PD_smp3','bonds.dump'), ...
  'bondTableFile', fullfile(root,'PD_smp3','bond.table'), ...
  'has_b_in_dump', true);   % this one does not need bond.table for N,b

Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'PD_smp4'), ...
  'bondsDumpFile', fullfile(root,'PD_smp4','bonds.dump'), ...
  'bondTableFile', fullfile(root,'PD_smp4','bond.table'), ...
  'has_b_in_dump', true);   % this one does not need bond.table for N,b

% Samples(end+1) = struct( ...
%   'dataDir',       fullfile(root,'PD_smp4'), ...
%   'bondsDumpFile', fullfile(root,'PD_smp4','bonds.dump'), ...
%   'bondTableFile', fullfile(root,'PD_smp4','bond.table'), ...
%   'has_b_in_dump', true);   % this one does not need bond.table for N,b

% Add more:
% Samples(end+1) = struct('dataDir', '...', 'bondsDumpFile','...\bonds.dump', ...
%                         'bondTableFile','...\bond.table', 'has_b_in_dump', true);

% --- Plot look ---
plotMode   = 'hist';            % 'hist' or 'line'
avg_mode   = 'pooled';          % 'pooled' or 'per-sample-mean'

% --- Consistent ranges (leave [] to auto from all samples) ---
r_range       = [];             % e.g., [0, 250]
N_range       = [];             % e.g., [1, 120] (integers)
lambda_range  = [];             % e.g., [0, 2]

% --- Bin counts (used for PDFs when 'hist' or 'line') ---
nbins_r       = 50;
nbins_lambda  = 50;             % continuous
% N uses integer-centered bins; edges built from N_range or data

% --- Misc ---
show_per_sample_overlays = false;
col_main = [86 184 112]/255;
color_line = [22 114 121]/255;

%% ===================== LOAD ALL SAMPLES (PASS 1) =====================
ns = numel(Samples);
assert(ns>=1, 'Please add at least one Sample in the "Samples" struct array.');

All_r = []; All_N = []; All_b = [];
PerSample = struct('r',[],'N',[],'b',[],'bond_type',[]);
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
        % Already provided in dump (N_in_dump, b_in_dump)
        keep = ~isnan(N_in_dump) & ~isnan(b_in_dump);
        if ~all(keep)
            warning('Sample %d: %d / %d bonds missing N,b in dump; dropping those.', s, nnz(~keep), numel(keep));
        end
        r   = r_d(keep);
        N   = N_in_dump(keep);
        b   = b_in_dump(keep);
        btp = bond_type(keep);
    else
        % Align with bond.table by pair keys
        [N_tab, b_tab] = read_bond_table_map(S.bondTableFile);
        [r, N, b, btp] = match_dump_with_table(i_d, j_d, r_d, bond_type, N_tab, b_tab);
    end

    lambda0 = r ./ (N .* b);

    % store per-sample
    PerSample(s).r = r;
    PerSample(s).N = N;
    PerSample(s).b = b;
    PerSample(s).lambda0 = lambda0;
    PerSample(s).bond_type = btp;

    All_r = [All_r; r]; %#ok<AGROW>
    All_N = [All_N; N]; %#ok<AGROW>
    All_b = [All_b; b]; %#ok<AGROW>

    fprintf('Sample %d: kept %d bonds. r[min,mean,max]=[%.3g, %.3g, %.3g], N[min,max]=[%g,%g], lambda[min,max]=[%.3g, %.3g]\n', ...
        s, numel(r), min(r), mean(r), max(r), min(N), max(N), min(lambda0), max(lambda0));
end

%% ===================== BUILD COMMON BINS =====================
% r (continuous)
if isempty(r_range)
    r_range = [min(All_r), max(All_r)];
end
edges_r = linspace(r_range(1), r_range(2), nbins_r+1);
xc_r    = 0.5*(edges_r(1:end-1) + edges_r(2:end));   % centers

% N (discrete integer-centered)
if isempty(N_range)
    N_range = [floor(min(All_N)), ceil(max(All_N))];
end
edgesN = (N_range(1)-0.5):(N_range(2)+0.5);  % integer-centered
xc_N   = N_range(1):N_range(2);

% lambda (continuous)
All_lambda = concatenate_field(PerSample,'lambda0');
if isempty(lambda_range)
    lambda_range = [min(All_lambda), max(All_lambda)];
end
edges_l = linspace(lambda_range(1), lambda_range(2), nbins_lambda+1);
xc_l    = 0.5*(edges_l(1:end-1) + edges_l(2:end));   % centers

%% ===================== HISTOGRAMS / PDFs =====================
switch avg_mode
    case 'pooled'
        r_pool = concatenate_field(PerSample,'r');
        N_pool = concatenate_field(PerSample,'N');
        l_pool = concatenate_field(PerSample,'lambda0');

        % --- per-bin probability (sum to 1), independent of bin width ---
        [~, prob_r] = binned_prob_per_bin_given_edges(r_pool, edges_r);
        [~, pmf_N ] = binned_pmf_discrete_given_edges(N_pool, edgesN);
        [~, prob_l] = binned_prob_per_bin_given_edges(l_pool, edges_l);

        % what we will plot
        y_r_mean = prob_r;
        y_N_mean = pmf_N;
        y_l_mean = prob_l;

    case 'per-sample-mean'
        Y_r = zeros(ns, numel(xc_r));
        Y_N = zeros(ns, numel(xc_N));
        Y_l = zeros(ns, numel(xc_l));
        for s = 1:ns
            [~, pr] = binned_prob_per_bin_given_edges(PerSample(s).r, edges_r);
            [~, pn] = binned_pmf_discrete_given_edges(PerSample(s).N, edgesN);
            [~, pl] = binned_prob_per_bin_given_edges(PerSample(s).lambda0, edges_l);
            Y_r(s,:) = pr;
            Y_N(s,:) = pn;
            Y_l(s,:) = pl;
        end
        y_r_mean = mean(Y_r,1);
        y_N_mean = mean(Y_N,1);
        y_l_mean = mean(Y_l,1);

    otherwise
        error('Unknown avg_mode: %s', avg_mode);
end


%% ===================== PLOTS =====================
figure('Name','End-to-end length r (averaged)'); hold on;
set(gca,'FontName','Times New Roman','FontSize',15);
if strcmp(plotMode,'hist')
    bar(xc_r, y_r_mean, 'FaceColor',col_main, 'EdgeColor',color_line, 'BarWidth',1);
else
    plot(xc_r, y_r_mean, 'LineWidth',1.8, 'Color', color_line);
end
xlabel('End-to-end length r','FontName','Times New Roman','FontSize',15);
ylabel('Probability per bin','FontName','Times New Roman','FontSize',15);
title(sprintf('r distribution (%s)', avg_mode),'FontName','Times New Roman','FontSize',15);
box on; grid off;

% ===================== N =====================
figure('Name','Kuhn segments N (averaged)'); hold on;
set(gca,'FontName','Times New Roman','FontSize',15);
if strcmp(plotMode,'hist')
    bar(xc_N, y_N_mean, 'FaceColor',col_main, 'EdgeColor',color_line, 'BarWidth',1);
else
    stem(xc_N, y_N_mean, 'LineWidth',1.6, 'Color', color_line, 'Marker','none');
end
xlabel('Kuhn segments N per bond','FontName','Times New Roman','FontSize',15);
ylabel('Probability per bin','FontName','Times New Roman','FontSize',15);
title(sprintf('N distribution (%s)', avg_mode),'FontName','Times New Roman','FontSize',15);
box on; grid off;

% ===================== lambda0 =====================
figure('Name','Pre-stretch \lambda_0 (averaged)'); hold on;
set(gca,'FontName','Times New Roman','FontSize',15);
if strcmp(plotMode,'hist')
    bar(xc_l, y_l_mean, 'FaceColor',col_main, 'EdgeColor',color_line, 'BarWidth',1);
else
    plot(xc_l, y_l_mean, 'LineWidth',1.8, 'Color', color_line);
end
xlabel('\lambda_0 = r / (N \cdot b)','FontName','Times New Roman','FontSize',15);
ylabel('Probability per bin','FontName','Times New Roman','FontSize',15);
title(sprintf('Pre-stretch distribution (%s)', avg_mode),'FontName','Times New Roman','FontSize',15);
box on; grid off;

% ===================== Overlays =====================
if show_per_sample_overlays && strcmp(avg_mode,'per-sample-mean')
    c = [0.2 0.2 0.2];
    for s = 1:ns
        [~, pr] = binned_prob_per_bin_given_edges(PerSample(s).r, edges_r);
        figure(findobj('Name','End-to-end length r (averaged)')); hold on;
        plot(xc_r, pr, ':', 'Color', c, 'LineWidth',1);

        [~, pn] = binned_pmf_discrete_given_edges(PerSample(s).N, edgesN);
        figure(findobj('Name','Kuhn segments N (averaged)')); hold on;
        plot(xc_N, pn, ':', 'Color', c, 'LineWidth',1);

        [~, pl] = binned_prob_per_bin_given_edges(PerSample(s).lambda0, edges_l);
        figure(findobj('Name','Pre-stretch \lambda_0 (averaged)')); hold on;
        plot(xc_l, pl, ':', 'Color', c, 'LineWidth',1);
    end
end



%% ===================== TEXT SUMMARY =====================
fprintf('\n=== GLOBAL SUMMARY (using common ranges) ===\n');
fprintf('r range used:       [%.3g, %.3g]\n', r_range(1), r_range(2));
fprintf('N range used:       [%d, %d]\n', N_range(1), N_range(2));
fprintf('lambda0 range used: [%.3g, %.3g]\n', lambda_range(1), lambda_range(2));
fprintf('Averaging mode: %s\n', avg_mode);
fprintf('sum(prob_r) = %.6f, sum(prob_l) = %.6f\n', sum(y_r_mean), sum(y_l_mean));