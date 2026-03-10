% Plot_Kuhn_Distributions_ExpFit.m
% MATLAB R2016a compatible (script, no functions)
%
% Goal:
  'bondTableFile', {}, ...
  'bondsDumpFile', {}, ...
  'label'        , {} );

root = 'E:\PhD\My Research\Polydisperse_fracture\PAPER';

Samples(end+1) = struct( ...
  'dataDir',       fullfile(root,'PD_smp2_notched'), ...
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


figure %Work of fracture
scatter([0 0 0 0],[6.41 6.28 6.35 6.42]*0.1,'o','LineWidth',2,'MarkerEdgeColor',[203 54 83]/255,'MarkerFaceColor',[221 132 146]/255)
hold on
scatter([5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49],[5.77 5.98 5.81 5.61 5.42 5.51 5.58 4.96 4.82 4.91 4.74 5.77 4.78 4.74 5.58 4.68 4.78 4.62 4.61 4.56 4.48 4.45]*0.1,'s','LineWidth',1,'MarkerEdgeColor',[22 114 121]/255,'MarkerFaceColor',[37 193 205]/255)
box on

figure %Toughness normalized
scatter((1/1.6)*[0 0 0 0],[4.8e4 4.8e4 5.1e4 5.02e4]/(4100*8.01),'o','LineWidth',1,'MarkerEdgeColor',[203 54 83]/255,'MarkerFaceColor',[221 132 146]/255)
hold on
scatter((1/1.6)*[5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49],[7.49e4 7.91e4 8.11e4 9.31e4 9.25e4 9.48e4 9.01e4 1.08e5 1.24e5 1.22e5 1.18e5 1.67e5 1.78e5 1.66e5 1.91e5 1.84e5 1.76e5 1.77e5 1.81e5 1.89e5 1.94e5 1.90e5]/(4100*12.2),'s','LineWidth',1,'MarkerEdgeColor',[22 114 121]/255,'MarkerFaceColor',[37 193 205]/255)
box on
L_dat = (1/1.6)*[0 0 0 0 5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49];
G_md = [4.8e4 4.8e4 5.1e4 5.02e4]/(4100*8.01);
G_pd = [7.49e4 7.91e4 8.11e4 9.31e4 9.25e4 9.48e4 9.01e4 1.08e5 1.24e5 1.22e5 1.18e5 1.67e5 1.78e5 1.66e5 1.91e5 1.84e5 1.76e5 1.77e5 1.81e5 1.89e5 1.94e5 1.90e5]/(4100*12.2);

G_dat = cat(2,G_md,G_pd);

AA = fit(L_dat',G_dat','poly1');
plot(AA,'k--')

figure %Fractocohesive length

l_fc_md = [4.8e4 4.8e4 5.1e4 5.02e4]/(4100)./([6.41 6.28 6.35 6.42]*0.1);
l_fc_pd = [7.49e4 7.91e4 8.11e4 9.31e4 9.25e4 9.48e4 9.01e4 1.08e5 1.24e5 1.22e5 1.18e5 1.67e5 1.78e5 1.66e5 1.91e5 1.84e5 1.76e5 1.77e5 1.81e5 1.89e5 1.94e5 1.90e5]/(4100)./([5.77 5.98 5.81 5.61 5.42 5.51 5.58 4.96 4.82 4.91 4.74 5.77 4.78 4.74 5.58 4.68 4.78 4.62 4.61 4.56 4.48 4.45]*0.1);

scatter((1/1.6)*[0 0 0 0],l_fc_md/23,'o','LineWidth',1,'MarkerEdgeColor',[203 54 83]/255,'MarkerFaceColor',[221 132 146]/255)
hold on
scatter((1/1.6)*[5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49],l_fc_pd/23,'s','LineWidth',1,'MarkerEdgeColor',[22 114 121]/255,'MarkerFaceColor',[37 193 205]/255)
l_fc_dat =  cat(2,l_fc_md,l_fc_pd)/23;

L_dat = (1/1.6)*[0 0 0 0 5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49];
H = fit(L_dat',l_fc_dat','poly1')
plot(H,'k--')
set(gca,'FontName','Times new roman')
set(gca,'FontSize',15)
box on


figure %standard deviation
L_dat = (1/1.6)*[0 0 0 0 5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49];
Sigma_dat = [0 0 0 0 4.87 4.87 4.87 9.45 9.45 9.45 9.45 14.19 14.19 14.19 14.19 20.45 20.45 20.45 20.45 25.19 25.19 25.19 25.19 30.42 30.42 30.42];
plot(L_dat,Sigma_dat);
ylabel('\sigma(L)')
xlabel('$\mathcal{L}/b$','Interpreter','latex')


figure %std divided by N0

N0_dat = [45 45 45 45 30 30 30 30 30 30 30 30 30 30 30 50 50 50 50 50 50 50 50 50 50 50];
plot(L_dat,Sigma_dat./N0_dat)
ylabel('\sigma(L)/N_0')
xlabel('$\mathcal{L}/b$','Interpreter','latex')

figure %std divided by mean
L_dat = (1/1.6)*[0 0 0 0 5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49];
Sigma_dat = [0 0 0 0 4.87 4.87 4.87 9.45 9.45 9.45 9.45 14.19 14.19 14.19 14.19 20.45 20.45 20.45 20.45 25.19 25.19 25.19 25.19 30.42 30.42 30.42];
Mean_dat = [45 45 45 45 30 30 30 40 40 40 40 45 45 45 45 70 70 70 70 75 75 75 75 80 80 80];

plot(L_dat,Sigma_dat./Mean_dat);
ylabel('\sigma(L)/<N>')
xlabel('$\mathcal{L}/b$','Interpreter','latex')


figure %C vs L

C = [-1.191 -1.188 -1.192 -1.208 -1.168 -1.175 -1.171 -1.185 -1.191 -1.202 -1.195 -1.091 -1.079 -1.095 -1.088 -1.002 -0.987 -0.994 -1.005 -1.004 -1.002 -0.991 -1.004 -1.376 -1.369 -1.381];
scatter(L_dat,exp(C))
ylabel('C')
xlabel('$\mathcal{L}/b$','Interpreter','latex')


x = [1 1.5 2 2.5 1 1.5 2 2.5 1 1.5 2 2.5 1 1.5 2 2.5 1 1.5 2 2.5];
y = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5];
c = [24.49 26.41 32.71 41.23 26.90 29.84 35.32 50.61 28.88 32.94 51.07 70.09 28.44 32.78 49.56 63.86 27.11 30.85 41.28 51.61];   % color-coded quantity

figure; hold on; box on;  %BD double sweep

scatter(x, y, 60, c, 'filled', ...
    'MarkerEdgeColor','k','LineWidth',0.75);

n = 256;
blue   = [0.0 0.2 0.8];
orange = [1.0 0.5 0.0];
red    = [0.8 0.0 0.0];

cmap = [ ...
    [linspace(blue(1),orange(1),n/2)' ...
     linspace(blue(2),orange(2),n/2)' ...
     linspace(blue(3),orange(3),n/2)']; ...
    [linspace(orange(1),red(1),n/2)' ...
     linspace(orange(2),red(2),n/2)' ...
     linspace(orange(3),red(3),n/2)'] ...
];

colormap(cmap);
colorbar;
caxis([min(c) max(c)]);

xlabel('x');
ylabel('y');
axis([0 4 0 6])

% --- Interpolation grid ---
xi = linspace(min(x), max(x), 200);
yi = linspace(min(y), max(y), 200);
[XI, YI] = meshgrid(xi, yi);

% --- Interpolate scattered c(x,y) ---
CI = griddata(x, y, c, XI, YI, 'natural');

% --- Contour levels ---
levels = linspace(min(c), max(c), 10);

% --- Overlay contour LINES ONLY ---
contour(XI, YI, CI, levels, ...
    'LineStyle','--', ...
    'LineColor','k', ...
    'LineWidth',1.0);


figure %Notch insensitivity

notch_lengths_pd = [32 32 32 80 80 80 160 160 160 320 320 320 640 640 640 870 870 870];
Work_fractures_pd = [1.172 1.245 1.194 1.231 1.201 1.188 0.924 0.899 0.914 0.688 0.697 0.694 0.567 0.540 0.532 0.456 0.481 0.449];

notch_lengths_md = [32 32 32 80 80 80 160 160 160 320 320 320 640 640 640 870 870 870];
Work_fractures_md = [2.151 2.117 2.091 1.910 1.887 1.814 1.325 1.277 1.394 0.732 0.797 0.794 0.708 0.680 0.702 0.681 0.649 0.699];


scatter(notch_lengths_pd/23,Work_fractures_pd,50,'s','LineWidth',1,'MarkerEdgeColor',[22 114 121]/255,'MarkerFaceColor',[37 193 205]/255)
axis([0 40 0 3])
hold on
scatter(notch_lengths_md/23,Work_fractures_md,50,'o','LineWidth',1,'MarkerEdgeColor',[203 54 83]/255,'MarkerFaceColor',[221 132 146]/255)

plot([100 100]/23,[0.5 1.5], 'LineWidth', 2, 'Color', [22 114 121]/255, 'LineStyle', ':');
plot([24 24]/23,[1.8 2.8],  'LineWidth', 2, 'Color', [203 54 88]/255, 'LineStyle', ':');


box on


figure %exp fit lengthscale

R0_md = [26.46 27.78 26.08 26.98];
R0_pd = [50.27 51.49 50.99 61.38 60.97 61.05 61.22 70.21 69.90 70.81 70.44 78.03 79.11 78.67 78.20 90.91 91.60 90.93 92.01 109.68 112.49 109.99];
R0_dat = cat(2,R0_md,R0_pd);
L_dat = (1/1.6)*[0 0 0 0 5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49];

scatter((1/1.6)*[0 0 0 0],(1/23)*R0_md,'o','LineWidth',2,'MarkerEdgeColor',[203 54 83]/255,'MarkerFaceColor',[221 132 146]/255)
hold on
scatter((1/1.6)*[5.12 5.12 5.12 10.21 10.21 10.21 10.21 16.46 16.46 16.46 16.46 20.52 20.52 20.52 20.52 25.41 25.41 25.41 25.41 30.49 30.49 30.49],(1/23)*R0_pd,'s','LineWidth',1,'MarkerEdgeColor',[22 114 121]/255,'MarkerFaceColor',[37 193 205]/255)
box on
F = fit(L_dat',(R0_dat/23)','poly1','lower',[])
plot(F,'k--')
% axis([0 20 30 100])

l_fc_md = (1/23)*[4.8e4 4.8e4 5.1e4 5.02e4]/(4100)./([6.41 6.28 6.35 6.42]*0.1);
l_fc_pd = (1/23)*[7.49e4 7.91e4 8.11e4 9.31e4 9.25e4 9.48e4 9.01e4 1.08e5 1.24e5 1.22e5 1.18e5 1.67e5 1.78e5 1.66e5 1.91e5 1.84e5 1.76e5 1.77e5 1.81e5 1.89e5 1.94e5 1.90e5]/(4100)./([5.77 5.98 5.81 5.61 5.42 5.51 5.58 4.96 4.82 4.91 4.74 5.77 4.78 4.74 5.58 4.68 4.78 4.62 4.61 4.56 4.48 4.45]*0.1);

figure
scatter((1/23)*R0_md,l_fc_md,'o','LineWidth',2,'MarkerEdgeColor',[203 54 83]/255,'MarkerFaceColor',[221 132 146]/255);
hold on
scatter((1/23)*R0_pd,l_fc_pd,'s','LineWidth',1,'MarkerEdgeColor',[22 114 121]/255,'MarkerFaceColor',[37 193 205]/255);
axis([0.5 5.5 0.5 5.5])
plot([0 120],[0 120],'k--','LineWidth',1.5)

l_fc_dat = cat(2,l_fc_md,l_fc_pd);
F = fit((1/23)*R0_dat',l_fc_dat','poly1')

plot(F,'k--')

xlabel('R_0')
ylabel('l_f_c')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',15)
set(gca,'FontSize',13);





