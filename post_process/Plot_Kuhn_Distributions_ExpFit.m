% Plot_Kuhn_Distributions_ExpFit.m
% MATLAB R2016a compatible (script, no functions)
%
% Goal:
%   - Read N_kuhn from bond.table for multiple samples
%   - Plot P(N) vs N for all samples on one figure
%   - Fit exponential tail P(N) ~ exp(-N/beta) for each sample
%     (beta = ? for monodisperse)
%   - ALSO: read bonds.dump (t=0) and compute mean pre-stretch
%           lambda0 = r/(N*b) for each sample

clear; clc; close all;

%% ===================== USER KNOBS =====================

% --- Samples list ---
% Each sample: directory + bond.table + bonds.dump path
Samples = struct( ...
  'dataDir'      , {}, ...
  'bondTableFile', {}, ...
  'bondsDumpFile', {}, ...
  'label'        , {} );

root = 'E:\PhD\My Research\Polydisperse_fracture\PAPER';

Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'PD_smp2_notched'), ...
  'bondTableFile', fullfile(root,'PD_smp2_notched','bondC.table'), ...
  'bondsDumpFile', fullfile(root,'PD_smp2_notched','bonds.dump'), ...
  'label',         'PD\_smp2');

Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'PD_widerdist_smp1_notched'), ...
  'bondTableFile', fullfile(root,'PD_widerdist_smp1_notched','bondC.table'), ...
  'bondsDumpFile', fullfile(root,'PD_widerdist_smp1_notched','bonds.dump'), ...
  'label',         'PD\_widest');

Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'PD_middist_smp1_notched'), ...
  'bondTableFile', fullfile(root,'PD_middist_smp1_notched','bondC.table'), ...
  'bondsDumpFile', fullfile(root,'PD_middist_smp1_notched','bonds.dump'), ...
  'label',         'PD\_wide');

Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'MD_smp2_notched'), ...
  'bondTableFile', fullfile(root,'MD_smp2_notched','bondC.table'), ...
  'bondsDumpFile', fullfile(root,'MD_smp2_notched','bonds.dump'), ...
  'label',         'MD\_smp1');

% --- Fitting options ---
% Minimum probability to include a bin in the exponential fit
p_min_fit = 1e-4;

%% ===================== LOAD N_kuhn + PRE-STRETCH =====================
ns  = numel(Samples);
assert(ns >= 1, 'Add at least one sample in the Samples struct.');

N_all_samples  = cell(ns,1);
All_N_global   = [];
lambda0_mean   = nan(ns,1);   % mean pre-stretch per sample

for s = 1:ns
    bt_file   = Samples(s).bondTableFile;
    bdump_file= Samples(s).bondsDumpFile;

    assert(exist(bt_file,'file')   == 2, 'Missing bond.table: %s', bt_file);
    assert(exist(bdump_file,'file')== 2, 'Missing bonds.dump: %s', bdump_file);

    fprintf('\n=== Sample %d/%d ===\n', s, ns);
    fprintf('  bond.table : %s\n', bt_file);
    fprintf('  bonds.dump : %s\n', bdump_file);

    % ---- parse bond.table ----
    % Expected numeric columns:
    %   1 index | 2 i | 3 j | 4 N_kuhn | 5 b_kuhn
    fid = fopen(bt_file,'r');
    if fid < 0
        error('Could not open %s', bt_file);
    end

    N_vec = [];

    tline = fgetl(fid);
    while ischar(tline)
        str = strtrim(tline);
        % skip empty or commented lines
        if ~isempty(str) && (str(1) ~= '#')
            % require the line to start with digit or sign
            if ~isempty(regexp(str,'^[\-\+]?\d','once'))
                toks = regexp(str,'\s+','split');
                vals = nan(1,5);
                nt   = numel(toks);
                m    = min(nt,5);
                j    = 1;
                while j <= m
                    v = str2double(toks{j});
                    if ~isnan(v)
                        vals(j) = v;
                    end
                    j = j + 1;
                end
                % keep if N_kuhn (4th col) is valid
                if ~isnan(vals(4))
                    N_vec = [N_vec; vals(4)]; %#ok<AGROW>
                end
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);

    if isempty(N_vec)
        warning('Sample %d: no N_kuhn found in bond.table, skipping P(N) for this sample.', s);
    else
        fprintf('  From bond.table: %d bonds. N[min, mean, max] = [%.3g, %.3g, %.3g]\n', ...
            numel(N_vec), min(N_vec), mean(N_vec), max(N_vec));
    end

    N_all_samples{s} = N_vec;
    All_N_global     = [All_N_global; N_vec]; %#ok<AGROW>

    % ---- parse bonds.dump at t=0 for pre-stretch ----
    % Expected helper: [bond_type, i_d, j_d, r_d, f_d, N_in_dump, b_in_dump]
    [bond_type, i_d, j_d, r_d, f_d, N_in_dump, b_in_dump] = read_bonds_dump_t0(bdump_file); %#ok<NASGU>

    keep = isfinite(r_d) & isfinite(N_in_dump) & isfinite(b_in_dump) & ...
           (N_in_dump > 0) & (b_in_dump > 0);

    if ~any(keep)
        warning('Sample %d: no valid bonds for pre-stretch in bonds.dump.', s);
        lambda0_mean(s) = NaN;
    else
        r0  = r_d(keep);
        N0  = N_in_dump(keep);
        b0  = b_in_dump(keep);
        lambda0_vec = r0 ./ (N0 .* b0);
        lambda0_mean(s) = mean(lambda0_vec);

        fprintf('  From bonds.dump t=0: mean pre-stretch lambda0 = %.4g\n', lambda0_mean(s));
    end
end

if isempty(All_N_global)
    error('No N_kuhn data found in any sample (bond.table).');
end

%% ===================== BUILD GLOBAL N GRID =====================
% N is integer -> use all distinct values observed
N_unique = unique(All_N_global);
nN       = numel(N_unique);

%% ===================== P(N) AND EXPONENTIAL FIT =====================
P_all  = zeros(ns, nN);   % P(N) per sample
beta   = nan(ns,1);       % exponential scale beta per sample

for s = 1:ns
    N_vec = N_all_samples{s};
    if isempty(N_vec)
        continue;
    end

    Nsamp = numel(N_vec);
    P_s   = zeros(1,nN);

    % empirical P(N)
    for k = 1:nN
        P_s(k) = sum(N_vec == N_unique(k)) / Nsamp;
    end

    P_all(s,:) = P_s;

    % detect monodisperse
    Nu = unique(N_vec);
    if numel(Nu) == 1
        beta(s) = Inf;
        fprintf('Sample %d (%s): monodisperse N = %.3g, beta = Inf\n', ...
            s, Samples(s).label, Nu);
        continue;
    end

    % -------- exponential fit: P(N) ~ exp(-N / beta) --------
    % Use bins with P > p_min_fit
    idx_fit = (P_s > p_min_fit);
    N_fit   = N_unique(idx_fit);
    P_fit   = P_s(idx_fit);

    if numel(N_fit) < 2
        warning('Sample %d: not enough points above p_min_fit for a robust fit.', s);
        beta(s) = NaN;
        continue;
    end

    y = log(P_fit(:));          % ln P
    X = [ones(numel(N_fit),1), N_fit(:)];  % [1  N]
    coeff = X \ y;              % least squares
    slope = coeff(2);

    beta(s) = -1 / slope;       % since ln P = const - N/beta

    fprintf('Sample %d (%s): fitted beta = %.4g (using %d points)\n', ...
        s, Samples(s).label, beta(s), numel(N_fit));
end

%% ===================== PLOT: P(N) vs N (all samples) =====================
figsize = [1 1 4 3];
set(0,'DefaultFigureUnits','inches','DefaultFigurePosition',figsize);
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',15);

colors = lines(ns);

figure('Color','w'); hold on; box on;
for s = 1:ns
    N_vec = N_all_samples{s};
    if isempty(N_vec)
        continue;
    end
    % only plot in [min(N_s), max(N_s)] for this sample
    Nmin = min(N_vec);
    Nmax = max(N_vec);
    mask = (N_unique >= Nmin) & (N_unique <= Nmax);

    plot(N_unique(mask), P_all(s,mask), '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 5, ...
        'Color', colors(s,:));
end

xlabel('Kuhn segment count N');
ylabel('P(N)');
title('Kuhn segment distributions P(N) for all samples');

leg_str = cell(ns,1);
for s = 1:ns
    if isinf(beta(s))
        leg_str{s} = sprintf('%s (\\beta = \\infty, \\lambda_0^{avg}=%.2f)', ...
                             Samples(s).label, lambda0_mean(s));
    else
        leg_str{s} = sprintf('%s (\\beta = %.2f, \\lambda_0^{avg}=%.2f)', ...
                             Samples(s).label, beta(s), lambda0_mean(s));
    end
end
legend(leg_str, 'Interpreter','tex', 'Location','northeastoutside');

%% (Optional) log-plot to visually check exponential tails
figure('Color','w'); hold on; box on;
for s = 1:ns
    N_vec = N_all_samples{s};
    if isempty(N_vec)
        continue;
    end
    Nmin = min(N_vec);
    Nmax = max(N_vec);
    mask = (N_unique >= Nmin) & (N_unique <= Nmax);

    semilogy(N_unique(mask), P_all(s,mask) + eps, '-o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 5, ...
        'Color', colors(s,:));
end
xlabel('Kuhn segment count N');
ylabel('P(N) (log scale)');
title('Kuhn segment distributions (semilogy)');
legend(leg_str, 'Interpreter','tex', 'Location','northeastoutside');


figure
scatter([0 0 0 0],[6.41 6.28 6.35 6.42]*0.1,'o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
hold on
scatter([16.46 16.46 16.46 16.46 20.52 20.52 20.52 30.49 30.49 30.49],[4.96 4.82 4.91 4.74 5.77 4.74 5.58 4.56 4.48 4.45]*0.1,'s','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[22 114 121]/255)

figure
scatter([0 0 0 0],[4.8e4 4.8e4 5.1e4 5.02e4]/(4100)./([6.41 6.28 6.35 6.42]*0.1),'o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
hold on
scatter([16.46 16.46 16.46 16.46 20.52 20.52 20.52 30.49 30.49 30.49],[7.9e4 8.5e4 8.8e4 9.1e4 1.67e5 1.66e5 1.91e5 2.08e5 2.04e5 2.02e5]/(4100)./([4.96 4.82 4.91 4.74 5.77 4.74 5.58 4.56 4.48 4.45]*0.1),'s','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[22 114 121]/255)




