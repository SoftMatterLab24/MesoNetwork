%% Post-processing: Stress–Stretch, W_frac, and beta from N_kuhn
% MATLAB R2016a compatible (script, no functions)
% Expected whitespace-delimited thermo columns:
%  1 Step | 2 v_ntotbonds | 3 c_virial[1] | 4 c_virial[2] | 5 c_virial[4] | 
%  6 PotEng | 7 KinEng | 8 v_kper | 9 c_fmax | 10 c_favg | 11 Lx | 12 Ly | 
% 13 c_virials[1] | 14 c_virials[2]

clear; clc; close all;

%% ---------------- User settings ----------------

% Choose which stress column to use for plotting sigma_yy:
%   'virial'  -> column 4 (c_virial[2])
%   'reduced' -> column 14 (c_virials[2])
stress_choice = 'virial';
% stress_choice = 'reduced';

% Optional: downsample points for plotting (set to 1 to keep all)
plot_stride = 1;

% Column index of N_kuhn in bond table (Bond.table / BondC.table)
N_kuhn_col = 4;   % adjust if your table format differs
nbins_beta = 25;  % number of bins for histogram in beta fit

% --------- File groups (fill these as needed) ---------
% THERMO FILES (as before)
mono_unnotched_files = { ...
%     'stress_md_smp1.txt', ...
    'stress_md_smp1.txt', ...
    'stress_md_smp1.txt', ...
    'stress_md_smp1.txt' , ...
    'stress_md_smp1.txt' , ... 
    };

mono_notched_files = {  ...
%     'stress_md_smp1_notched.txt', ...
    'stress_md_smp2_notched.txt', ...
    'stress_md_smp3_notched.txt', ...
    'stress_md_smp4_notched.txt',...
    'stress_md_smp5_notched.txt'
    };

poly_unnotched_files = { ...
    'stress_1e-6.txt', ...
    'stress_smp2.txt', ...
    'stress_smp3.txt', ...
    'stress_smp4.txt', ...
    'stress_pd_widerdist_smp1.txt', ...
    'stress_pd_widerdist_smp2.txt', ...
    'stress_pd_widerdist_smp3.txt', ...
    'stress_pd_middist_smp1.txt', ...
    'stress_pd_middist_smp2.txt', ...
    'stress_pd_middist_smp3.txt', ...
};

poly_notched_files = { ...
    'stress_1e-6_notched.txt', ...
    'stress_smp2_notched.txt', ...
    'stress_smp3_notched.txt', ...
    'stress_smp4_notched.txt', ...
    'stress_pd_widerdist_smp1_notched.txt', ...
    'stress_pd_widerdist_smp2_notched.txt', ...
    'stress_pd_widerdist_smp2_notched.txt', ...
    'stress_pd_middist_smp1_notched.txt', ...
    'stress_pd_middist_smp2_notched.txt', ...
    'stress_pd_middist_smp3_notched.txt', ...
};

% --------- NEW: directories that contain bondC.table for beta ---------
% Each directory must contain a file named exactly 'bondC.table'.
% Order must correspond to the UNNOTCHED thermo files used in the pairs.

mono_bond_dirs = { ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp1', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp3_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp4_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp5_notched', ...
};

poly_bond_dirs = { ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\sample_alpine_huge_erate_1e-6_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp3_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp4_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp1_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp3_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_middist_smp1_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_middist_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_middist_smp3_notched', ...
};



% --------- Flatten thermo groups into a single list + metadata ---------
files      = {};
isMono     = [];
isNotched  = [];
all_labels = {};

%%% keep track of index pairs for MD and PD
mono_un_idx    = [];
mono_notch_idx = [];
poly_un_idx    = [];
poly_notch_idx = [];

% mono, unnotched
for i = 1:numel(mono_unnotched_files)
    f = mono_unnotched_files{i};
    files{end+1} = f; %#ok<SAGROW>
    isMono(end+1)    = 1; %#ok<SAGROW>
    isNotched(end+1) = 0; %#ok<SAGROW>
    mono_un_idx(end+1) = numel(files); %#ok<SAGROW>
end

% mono, notched
for i = 1:numel(mono_notched_files)
    f = mono_notched_files{i};
    files{end+1} = f; %#ok<SAGROW>
    isMono(end+1)    = 1; %#ok<SAGROW>
    isNotched(end+1) = 1; %#ok<SAGROW>
    mono_notch_idx(end+1) = numel(files); %#ok<SAGROW>
end

% poly, unnotched
for i = 1:numel(poly_unnotched_files)
    f = poly_unnotched_files{i};
    files{end+1} = f; %#ok<SAGROW>
    isMono(end+1)    = 0; %#ok<SAGROW>
    isNotched(end+1) = 0; %#ok<SAGROW>
    poly_un_idx(end+1) = numel(files); %#ok<SAGROW>
end

% poly, notched
for i = 1:numel(poly_notched_files)
    f = poly_notched_files{i};
    files{end+1} = f; %#ok<SAGROW>
    isMono(end+1)    = 0; %#ok<SAGROW>
    isNotched(end+1) = 1; %#ok<SAGROW>
    poly_notch_idx(end+1) = numel(files); %#ok<SAGROW>
end

if isempty(files)
    error('No input thermo files specified in the four groups.');
end

% Build labels from filenames
for k = 1:numel(files)
    [~, base, ext] = fileparts(files{k});
    all_labels{k,1} = [base, ext];
end

nF = numel(files);

% Preallocate cell arrays for results
all_lambda_y = cell(nF,1);
all_sigma_yy = cell(nF,1);
all_KinEng   = cell(nF,1);  % still read, but not plotted
all_kper     = cell(nF,1);
all_fmax     = cell(nF,1);

%% --------------- Parse each thermo file -----------------
for k = 1:nF
    fname = files{k};
    fid = fopen(fname,'r');
    if fid < 0
        error('Cannot open file: %s', fname);
    end
    
    rows = [];
    
    % Read line-by-line to safely skip headers and non-numeric lines
    tline = fgetl(fid);
    while ischar(tline)
        s = strtrim(tline);
        if ~isempty(s)
            % Accept rows that start with a digit or sign
            if ~isempty(regexp(s, '^[\-\+]?\d', 'once'))
                toks = regexp(s, '\s+', 'split');
                vals = nan(1,14);
                nt = numel(toks);
                m  = min(nt,14);
                j = 1;
                while j <= m
                    v = str2double(toks{j});
                    if ~isnan(v)
                        vals(j) = v;
                    end
                    j = j + 1;
                end
                if ~isnan(vals(1)) && ~isnan(vals(12))
                    rows = [rows; vals]; %#ok<AGROW>
                end
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    
    if isempty(rows)
        error('No numeric thermo rows found in %s', fname);
    end
    
    % Map columns
    Step   = rows(:,1);
    cvir2  = rows(:,4);   % yy (raw virial)
    PotEng = rows(:,6); %#ok<NASGU>
    KinEng = rows(:,7);
    kper   = rows(:,8);
    fmax   = rows(:,9);
    Ly     = rows(:,12);
    cvirs2 = rows(:,14);  % yy (reduced/alt)
    
    % Remove rows with NaNs in key columns
    good = isfinite(Step) & isfinite(Ly) & isfinite(KinEng) & ...
           isfinite(kper) & isfinite(fmax) & (isfinite(cvir2) | isfinite(cvirs2));
    Step   = Step(good);
    cvir2  = cvir2(good);
    cvirs2 = cvirs2(good);
    KinEng = KinEng(good);
    kper   = kper(good);
    fmax   = fmax(good);
    Ly     = Ly(good);
    
    % Deduplicate by Step (keep last)
    [~, idxLast] = unique(Step, 'last');
    keep = false(size(Step));
    keep(idxLast) = true;
    
    Step   = Step(keep);
    cvir2  = cvir2(keep);
    cvirs2 = cvirs2(keep);
    KinEng = KinEng(keep);
    kper   = kper(keep);
    fmax   = fmax(keep);
    Ly     = Ly(keep);
    
    % Sort by Step increasing
    [Step, ord] = sort(Step);
    cvir2  = cvir2(ord);
    cvirs2 = cvirs2(ord);
    KinEng = KinEng(ord);
    kper   = kper(ord);
    fmax   = fmax(ord);
    Ly     = Ly(ord);
    
    % Compute stretch in Y (loading direction)
    if isempty(Ly) || Ly(1) == 0
        error('Ly(1) is empty or zero in %s; cannot compute lambda_y.', fname);
    end
    lambda_y = Ly ./ Ly(1);
    
    % Choose stress column
    switch lower(stress_choice)
        case 'virial'
            sigma_yy = cvir2;
        case 'reduced'
            sigma_yy = cvirs2;
        otherwise
            error('Unknown stress_choice: %s', stress_choice);
    end
    
    % Downsample for plotting if requested
    if plot_stride > 1
        idx = 1:plot_stride:numel(lambda_y);
        lambda_y = lambda_y(idx);
        sigma_yy = sigma_yy(idx);
        KinEng   = KinEng(idx);
        kper     = kper(idx);
        fmax     = fmax(idx);
    end
    
    all_lambda_y{k} = lambda_y;
    all_sigma_yy{k} = sigma_yy;
    all_KinEng{k}   = KinEng;
    all_kper{k}     = kper;
    all_fmax{k}     = fmax;
end

%% --------------- Fracture energy (per sample, paired) --------------
% For each sample pair (unnotched, notched):
%   1) From NOTCHED: find lambda_fail where sigma drops to 0.4*sigma_max
%   2) From UNNOTCHED: integrate -sigma_yy vs lambda_y from 1 to lambda_fail

lambda_fail_notched = nan(nF,1);

% --- Step 1: failure stretch from all notched curves (0.4*sigma_max) ---
for k = 1:nF
    if isNotched(k)
        lam = all_lambda_y{k};
        sig = -all_sigma_yy{k};  % positive tension
        
        if ~isempty(lam) && ~isempty(sig)
            n_local = min(numel(lam), numel(sig));
            lam = lam(1:n_local);
            sig = sig(1:n_local);
            
            [sig_max, imax] = max(sig);
            if ~isempty(imax) && isfinite(sig_max)
                sig_post = sig(imax:end);
                lam_post = lam(imax:end);
                
                thresh   = 0.4 * sig_max;
                idx_drop = find(sig_post <= thresh, 1, 'first');
                
                if ~isempty(idx_drop)
                    lambda_fail_notched(k) = lam_post(idx_drop);
                else
                    lambda_fail_notched(k) = lam(imax);
                end
            end
        end
    end
end

% --- Step 2: MD and PD pairs ---
nMonoPairs = min(numel(mono_un_idx),  numel(mono_notch_idx));
nPolyPairs = min(numel(poly_un_idx),  numel(poly_notch_idx));

% Sanity check (optional)
if numel(mono_bond_dirs) < nMonoPairs
    warning('mono_bond_dirs has fewer entries (%d) than mono pairs (%d).', ...
        numel(mono_bond_dirs), nMonoPairs);
end
if numel(poly_bond_dirs) < nPolyPairs
    warning('poly_bond_dirs has fewer entries (%d) than poly pairs (%d).', ...
        numel(poly_bond_dirs), nPolyPairs);
end



lambda_fail_mono = nan(nMonoPairs,1);
lambda_fail_poly = nan(nPolyPairs,1);
W_mono           = nan(nMonoPairs,1);
W_poly           = nan(nPolyPairs,1);

fprintf('\n=== Per-sample failure stretch and work of fracture ===\n');

% ----- MONODISPERSE pairs -----
if nMonoPairs > 0
    fprintf('\n--- Monodisperse samples ---\n');
    for j = 1:nMonoPairs
        k_un   = mono_un_idx(j);
        k_not  = mono_notch_idx(j);
        lam_un = all_lambda_y{k_un};
        sig_un = -all_sigma_yy{k_un};   % positive tension
        lf     = lambda_fail_notched(k_not);
        
        lambda_fail_mono(j) = lf;
        
        if ~isnan(lf)
            mask = (lam_un >= 1) & (lam_un <= lf);
            lam_clip = lam_un(mask);
            sig_clip = sig_un(mask);
            if numel(lam_clip) >= 2
                W_mono(j) = trapz(lam_clip, sig_clip);
            end
        end
        
        fprintf('  Pair %d: UN=%s | NOTCH=%s\n', j, ...
            all_labels{k_un}, all_labels{k_not});
        fprintf('      lambda_fail = %.4f,  W = %.4e\n', ...
            lambda_fail_mono(j), W_mono(j));
    end
    
    if any(isfinite(W_mono))
        fprintf('  --> Mean W_mono = %.4e\n', mean(W_mono(isfinite(W_mono))));
    end
end

% ----- POLYDISPERSE pairs -----
if nPolyPairs > 0
    fprintf('\n--- Polydisperse samples ---\n');
    for j = 1:nPolyPairs
        k_un   = poly_un_idx(j);
        k_not  = poly_notch_idx(j);
        lam_un = all_lambda_y{k_un};
        sig_un = -all_sigma_yy{k_un};   % positive tension
        lf     = lambda_fail_notched(k_not);
        
        lambda_fail_poly(j) = lf;
        
        if ~isnan(lf)
            mask = (lam_un >= 1) & (lam_un <= lf);
            lam_clip = lam_un(mask);
            sig_clip = sig_un(mask);
            if numel(lam_clip) >= 2
                W_poly(j) = trapz(lam_clip, sig_clip);
            end
        end
        
        fprintf('  Pair %d: UN=%s | NOTCH=%s\n', j, ...
            all_labels{k_un}, all_labels{k_not});
        fprintf('      lambda_fail = %.4f,  W = %.4e\n', ...
            lambda_fail_poly(j), W_poly(j));
    end
    
    if any(isfinite(W_poly))
        fprintf('  --> Mean W_poly = %.4e\n', mean(W_poly(isfinite(W_poly))));
    end
end

%% ------------------ Compute beta from bond tables ---------------------

%% ------------------ Compute beta from bond tables ---------------------

% Monodisperse: all chains have same N_kuhn -> beta = 0 by definition
beta_mono = zeros(nMonoPairs,1);

% Polydisperse: compute beta from bondC.table in poly_bond_dirs
beta_poly = nan(nPolyPairs,1);

fprintf('\n=== Beta from N_kuhn distributions (poly only) ===\n');

for j = 1:nPolyPairs
    if j > numel(poly_bond_dirs)
        warning('No poly_bond_dirs entry for pair %d.', j);
        continue;
    end
    
    btfile = fullfile(poly_bond_dirs{j}, 'bondC.table');
    fprintf('\n[POLY] Pair %d: reading %s\n', j, btfile);
    
    fid = fopen(btfile,'r');
    if fid < 0
        warning('Cannot open bond table: %s', btfile);
        continue;
    end
    
    Nvals = [];
    tline = fgetl(fid);
    while ischar(tline)
        s = strtrim(tline);
        if ~isempty(s)
            % numeric lines only
            if ~isempty(regexp(s,'^[\-\+]?\d','once'))
                toks = regexp(s,'\s+','split');
                nt   = numel(toks);
                if nt >= N_kuhn_col
                    vN = str2double(toks{N_kuhn_col});
                    if ~isnan(vN)
                        Nvals(end+1,1) = vN; %#ok<AGROW>
                    end
                end
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    
    Nvals = Nvals(Nvals > 0);
    fprintf('  Pair %d: %d N_kuhn values after >0 filter\n', j, numel(Nvals));
    
    if numel(Nvals) < 10
        warning('Not enough N_kuhn data in %s for a stable beta fit.', btfile);
        continue;
    end
    
    % Histogram + exponential fit: P(N) ~ exp(-N/beta)
    Nmin = min(Nvals);
    Nmax = max(Nvals);
    if Nmax <= Nmin
        warning('Nmax <= Nmin in %s; skipping.', btfile);
        continue;
    end
    
    edges   = linspace(Nmin, Nmax, nbins_beta+1);
    counts  = histcounts(Nvals, edges);
    centers = 0.5*(edges(1:end-1) + edges(2:end));
    P       = counts / sum(counts);
    mask    = (P > 0);
    xfit    = centers(mask);
    yfit    = P(mask);
    
    if numel(xfit) < 3
        warning('Too few non-zero bins for beta fit in %s.', btfile);
        continue;
    end
    
    logy = log(yfit);
    p = polyfit(xfit, logy, 1);   % log P ~ a + b N
    b = p(1);
    fprintf('  Fitted slope b = %.4e\n', b);
    
    if b >= 0
        warning('Fitted slope b >= 0 in %s; beta undefined.', btfile);
        continue;
    end
    
    beta_poly(j) = -1 / b;
    fprintf('  ==> Poly pair %d: beta = %.4f\n', j, beta_poly(j));
end



%% ------------------ Plotting: Stress–Stretch ---------------------
figsize = [1 1 3 3];
set(0,'DefaultFigureUnits','inches','DefaultFigurePosition',figsize);
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',15);
cMono = [203 54 88]/255;
cPoly = [22  114 121]/255;

% Monodisperse + Polydisperse in one figure (like before)
idx_mono = find(isMono ~= 0);
idx_poly = find(isMono == 0);

figure('Color','w'); hold on; box on;

% Monodisperse
if ~isempty(idx_mono)
    for j = 1:numel(idx_mono)
        k = idx_mono(j);
        if isNotched(k)
            ls = '--';
        else
            ls = '-';
        end
        plot(all_lambda_y{k}, -all_sigma_yy{k}, ls, ...
            'LineWidth', 1.5, 'Color', cMono);
    end
end

% Polydisperse
if ~isempty(idx_poly)
    for j = 1:numel(idx_poly)
        k = idx_poly(j);
        if isNotched(k)
            ls = '--';
        else
            ls = '-';
        end
        plot(all_lambda_y{k}, -all_sigma_yy{k}, ls, ...
            'LineWidth', 1.5, 'Color', cPoly);
    end
end

xlabel('\lambda_y','Interpreter','tex');
ylabel('\sigma_{yy}','Interpreter','tex');
title('\sigma_{yy} vs. \lambda_y');
legend(all_labels, 'Interpreter','none', 'Location','best');
axis([1 4.5 0 4]);

%% ------------------ Plotting: Notched-only stress–stretch --------
idx_notched = find(isNotched ~= 0);

if ~isempty(idx_notched)
    figure('Color','w'); hold on; box on;
    
    for j = 1:numel(idx_notched)
        k = idx_notched(j);
        if isMono(k)
            col = cMono;
        else
            col = cPoly;
        end
        plot(all_lambda_y{k}, -all_sigma_yy{k}, '-', ...
            'LineWidth', 2, 'Color', col);
    end
    
    xlabel('\lambda_y','Interpreter','tex');
    ylabel('\sigma_{yy}','Interpreter','tex');
    title('Notched samples only: \sigma_{yy} vs. \lambda_y');
    legend(all_labels(idx_notched), 'Interpreter','none', 'Location','best');
    axis([1 3 0 1.5]);
end

%% ------------------ Plotting: W_frac vs beta ---------------------
figure('Color','w'); hold on; box on;

% Mono scatter
maskM = isfinite(beta_mono) & isfinite(W_mono);
if any(maskM)
    plot(beta_mono(maskM), W_mono(maskM), 's', ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerFaceColor', cMono, 'MarkerEdgeColor', cMono);
end

% Poly scatter
maskP = isfinite(beta_poly) & isfinite(W_poly);
if any(maskP)
    plot(beta_poly(maskP), W_poly(maskP), 'o', ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerFaceColor', cPoly, 'MarkerEdgeColor', cPoly);
end

xlabel('\beta','Interpreter','tex');
ylabel('W_{\mathrm{frac}}','Interpreter','tex');
title('Work of fracture vs. polydispersity \beta');
leg = {};
if any(maskM), leg{end+1} = 'Monodisperse pairs'; end
if any(maskP), leg{end+1} = 'Polydisperse pairs'; end
if ~isempty(leg)
    legend(leg, 'Location','best');
end
