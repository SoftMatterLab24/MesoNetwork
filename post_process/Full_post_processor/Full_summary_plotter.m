% -------------------------------------------------------------------------
% SummaryPlotter_FromPairSummaries.m
%
% Re-plot global quantities from per-pair summaries stored in NOTCHED dirs.
%
% For each NOTCHED directory:
%   - Reads stress*.txt    -> notched stress–stretch (sigma_yy vs lambda_y)
%   - Reads *_pair_summary.txt -> beta, W_frac, Enorm, l_fc, group, label
%
% Global plots:
%   - Notched sigma_yy vs lambda_y
%   - W_frac vs beta
%   - Enorm vs beta
%   - l_fc vs beta
%
% Assumes pair summary files were written by the main post-processing script:
%   Pair label       : <string>
%   Group            : mono | poly | bimo
%   ...
%   beta             : <number>      [-]
%   W_frac           : <number>      [MPa]
%   Enorm            : <number>      [-]
%   l_fc             : <number>      [nm]
%
% Requirements:
%   - MATLAB R2016a compatible (script, no local functions, no "contains").
% -------------------------------------------------------------------------
clear; clc; close all;

%% ===================== USER INPUT: NOTCHED DIR LISTS ===================
% Each cell array below contains ONLY *NOTCHED* folders.
% The script finds the corresponding *_pair_summary.txt inside each folder.

% ----------------- MONODISPERSE notched dirs -----------------
monoNotchedDirs = {
   'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp2_notched';
};

% ----------------- POLYDISPERSE notched dirs -----------------
polyNotchedDirs = {
   'E:\PhD\My Research\Polydisperse_fracture\PAPER\sample_alpine_huge_erate_1e-6_notched';
   'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp1_notched';
};

% ----------------- BIMODAL notched dirs -----------------
bimoNotchedDirs = {
%    'E:\...\BD_smp1_notched';
};

% Thermo file name pattern inside each NOTCHED directory:
thermoPattern = 'stress*.txt';

% Choose which stress column to use for plotting sigma_yy:
%   'virial'  -> column 4  (c_virial[2])
%   'reduced' -> column 14 (c_virials[2])
stress_choice = 'virial';
% stress_choice = 'reduced';

% Optional downsample for plotting (set to 1 to keep all points)
plot_stride = 1;

%% ===================== BUILD LIST OF ALL NOTCHED DIRS ==================
notDirs   = {};
groups    = {};
dirLabels = {};  % will be overwritten by pair summary "Pair label" if present

% Mono
for i = 1:numel(monoNotchedDirs)
    notDirs{end+1,1}   = monoNotchedDirs{i};
    groups{end+1,1}    = 'mono';
    dirLabels{end+1,1} = sprintf('Mono_%d', i);
end

% Poly
for i = 1:numel(polyNotchedDirs)
    notDirs{end+1,1}   = polyNotchedDirs{i};
    groups{end+1,1}    = 'poly';
    dirLabels{end+1,1} = sprintf('Poly_%d', i);
end

% Bimo
for i = 1:numel(bimoNotchedDirs)
    notDirs{end+1,1}   = bimoNotchedDirs{i};
    groups{end+1,1}    = 'bimo';
    dirLabels{end+1,1} = sprintf('Bimo_%d', i);
end

nPairs = numel(notDirs);
if nPairs == 0
    error('No notched directories specified. Fill monoNotchedDirs / polyNotchedDirs / bimoNotchedDirs.');
end

%% ===================== STORAGE FOR DATA ================================
beta_all   = nan(nPairs,1);
Wfrac_all  = nan(nPairs,1);
Enorm_all  = nan(nPairs,1);
lfc_all    = nan(nPairs,1);
group_all  = groups;       % cell of strings
label_all  = dirLabels;    % will be updated from summary if available

all_lambda_not = cell(nPairs,1);
all_sigma_not  = cell(nPairs,1);

%% ===================== LOOP OVER NOTCHED DIRS ==========================
fprintf('\n=========== Summary Plotter: %d notched dirs ===========\n', nPairs);

for p = 1:nPairs
    notDir = notDirs{p};
    group_guess = groups{p};
    
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Pair %d / %d\n', p, nPairs);
    fprintf('  Notched dir (input) : %s\n', notDir);
    fprintf('  Group (from input)  : %s\n', group_guess);
    
    %% ---- 1) Find and read pair summary txt ----------------------------
    d_sum = dir(fullfile(notDir, '*_pair_summary.txt'));
    if isempty(d_sum)
        warning('  No *_pair_summary.txt found in %s. Skipping beta/W_frac/Enorm/l_fc.', notDir);
    else
        sumFile = fullfile(notDir, d_sum(1).name);
        fprintf('  Using summary file   : %s\n', sumFile);
        
        fidS = fopen(sumFile, 'r');
        if fidS == -1
            warning('  Could not open summary file. Skipping scalars.');
        else
            % Initialize local variables
            label_local = label_all{p};
            group_local = group_guess;
            beta_local  = NaN;
            Wfrac       = NaN;
            Enorm       = NaN;
            lfc         = NaN;
            
            line = fgetl(fidS);
            while ischar(line)
                s = strtrim(line);
                
                % Pair label
                if ~isempty(s) && ~isempty(strfind(s, 'Pair label'))
                    tok = regexp(s, '^Pair label\s*:\s*(.*)$', 'tokens', 'once');
                    if ~isempty(tok)
                        label_local = strtrim(tok{1});
                    end
                end
                
                % Group
                if ~isempty(s) && ~isempty(strfind(s, 'Group'))
                    tok = regexp(s, '^Group\s*:\s*(.*)$', 'tokens', 'once');
                    if ~isempty(tok)
                        group_local = strtrim(tok{1});
                    end
                end
                
                % beta
                if ~isempty(s) && strncmp(s, 'beta', 4)
                    tok = regexp(s, 'beta\s*:\s*([Ee0-9\+\-\.]+)', 'tokens', 'once');
                    if ~isempty(tok)
                        beta_local = str2double(tok{1});
                    end
                end
                
                % W_frac
                if ~isempty(s) && strncmp(s, 'W_frac', 6)
                    tok = regexp(s, 'W_frac\s*:\s*([Ee0-9\+\-\.]+)', 'tokens', 'once');
                    if ~isempty(tok)
                        Wfrac = str2double(tok{1});
                    end
                end
                
                % Enorm
                if ~isempty(s) && strncmp(s, 'Enorm', 5)
                    tok = regexp(s, 'Enorm\s*:\s*([Ee0-9\+\-\.]+)', 'tokens', 'once');
                    if ~isempty(tok)
                        Enorm = str2double(tok{1});
                    end
                end
                
                % l_fc
                if ~isempty(s) && strncmp(s, 'l_fc', 4)
                    tok = regexp(s, 'l_fc\s*:\s*([Ee0-9\+\-\.]+)', 'tokens', 'once');
                    if ~isempty(tok)
                        lfc = str2double(tok{1});
                    end
                end
                
                line = fgetl(fidS);
            end
            fclose(fidS);
            
            % Store back
            label_all{p} = label_local;
            group_all{p} = group_local;
            beta_all(p)  = beta_local;
            Wfrac_all(p) = Wfrac;
            Enorm_all(p) = Enorm;
            lfc_all(p)   = lfc;
            
            fprintf('  Parsed from summary: label=%s, group=%s, beta=%.4g, W_frac=%.4g, Enorm=%.4g, l_fc=%.4g\n', ...
                label_local, group_local, beta_local, Wfrac, Enorm, lfc);
        end
    end
    
    %% ---- 2) Read NOTCHED thermo for sigma_yy vs lambda_y --------------
    d_not = dir(fullfile(notDir, thermoPattern));
    if isempty(d_not)
        warning('  No thermo file matching "%s" in %s. Skipping stress-stretch.', ...
            thermoPattern, notDir);
        continue;
    end
    thermo_not = fullfile(notDir, d_not(1).name);
    fprintf('  Notched thermo file : %s\n', thermo_not);
    
    fidT = fopen(thermo_not, 'r');
    if fidT < 0
        warning('  Cannot open thermo file: %s. Skipping stress-stretch.', thermo_not);
        continue;
    end
    
    rows = [];
    tline = fgetl(fidT);
    while ischar(tline)
        s = strtrim(tline);
        if ~isempty(s)
            % numeric line starts with +/- digit
            if ~isempty(regexp(s, '^[\-\+]?\d', 'once'))
                toks = regexp(s, '\s+', 'split');
                vals = nan(1,14);
                nt   = numel(toks);
                m    = min(nt,14);
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
        tline = fgetl(fidT);
    end
    fclose(fidT);
    
    if isempty(rows)
        warning('  No numeric rows in thermo file for %s. Skipping stress-stretch.', notDir);
        continue;
    end
    
    Step_n   = rows(:,1);
    cvir2_n  = rows(:,4);
    cvirs2_n = rows(:,14);
    Ly_n     = rows(:,12);
    
    good = isfinite(Step_n) & isfinite(Ly_n) & ...
           (isfinite(cvir2_n) | isfinite(cvirs2_n));
    Step_n   = Step_n(good);
    cvir2_n  = cvir2_n(good);
    cvirs2_n = cvirs2_n(good);
    Ly_n     = Ly_n(good);
    
    [~, idxLast] = unique(Step_n, 'last');
    keep = false(size(Step_n));
    keep(idxLast) = true;
    
    Step_n   = Step_n(keep);
    cvir2_n  = cvir2_n(keep);
    cvirs2_n = cvirs2_n(keep);
    Ly_n     = Ly_n(keep);
    
    [Step_n, ord] = sort(Step_n);
    cvir2_n  = cvir2_n(ord);
    cvirs2_n = cvirs2_n(ord);
    Ly_n     = Ly_n(ord);
    
    if isempty(Ly_n) || Ly_n(1) == 0
        warning('  Ly(1) invalid for %s. Skipping stress-stretch.', notDir);
    else
        lambda_not = Ly_n ./ Ly_n(1);
        
        switch lower(stress_choice)
            case 'virial'
                sigma_not_raw = cvir2_n;
            case 'reduced'
                sigma_not_raw = cvirs2_n;
            otherwise
                error('Unknown stress_choice: %s', stress_choice);
        end
        sigma_not = -sigma_not_raw;   % tension positive
        
        if plot_stride > 1
            idx = 1:plot_stride:numel(lambda_not);
            lambda_not = lambda_not(idx);
            sigma_not  = sigma_not(idx);
        end
        
        all_lambda_not{p} = lambda_not;
        all_sigma_not{p}  = sigma_not;
    end
end

%% ===================== PLOTTING SETTINGS ===============================
figsize = [1 1 4 3];
set(0,'DefaultFigureUnits','inches','DefaultFigurePosition',figsize);
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',14);

cMono = [203 54 88]/255;
cPoly = [22 114 121]/255;
cBimo = [0.2 0.2 0.2];

%% ===================== PLOT: NOTCHED STRESS–STRETCH ====================
figure('Color','w'); hold on; box on;

for p = 1:nPairs
    lam_not = all_lambda_not{p};
    sig_not = all_sigma_not{p};
    if isempty(lam_not) || isempty(sig_not)
        continue;
    end
    
    switch lower(group_all{p})
        case 'mono'
            col = cMono;
        case 'poly'
            col = cPoly;
        case 'bimo'
            col = cBimo;
        otherwise
            col = [0 0 0];
    end
    
    plot(lam_not, sig_not, '-', 'LineWidth', 1.5, 'Color', col);
end

xlabel('\lambda_y','Interpreter','tex');
ylabel('\sigma_{yy}','Interpreter','tex');
title('\sigma_{yy} vs. \lambda_y (NOTCHED only)','Interpreter','tex');
box on;

%% ===================== PLOT: W_frac vs beta ============================
figure('Color','w'); hold on; box on;

for p = 1:nPairs
    if ~isfinite(beta_all(p)) || ~isfinite(Wfrac_all(p))
        continue;
    end
    
    switch lower(group_all{p})
        case 'mono'
            marker = 's';
            col    = cMono;
        case 'poly'
            marker = 'o';
            col    = cPoly;
        case 'bimo'
            marker = 'd';
            col    = cBimo;
        otherwise
            marker = '^';
            col    = [0 0 0];
    end
    
    plot(beta_all(p), Wfrac_all(p), marker, ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
end

xlabel('\beta','Interpreter','tex');
ylabel('W_{frac}','Interpreter','tex');   % [MPa]
title('Work of fracture vs. \beta','Interpreter','tex');
legend({'Mono','Poly','Bimodal'}, 'Location','best');
box on;

%% ===================== PLOT: Enorm vs beta =============================
figure('Color','w'); hold on; box on;

for p = 1:nPairs
    if ~isfinite(beta_all(p)) || ~isfinite(Enorm_all(p))
        continue;
    end
    
    switch lower(group_all{p})
        case 'mono'
            marker = 's';
            col    = cMono;
        case 'poly'
            marker = 'o';
            col    = cPoly;
        case 'bimo'
            marker = 'd';
            col    = cBimo;
        otherwise
            marker = '^';
            col    = [0 0 0];
    end
    
    plot(beta_all(p), Enorm_all(p), marker, ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
end

xlabel('\beta','Interpreter','tex');
ylabel('E_{norm}','Interpreter','tex');   % = G_diss / T0
title('Normalized rupture energy vs. \beta','Interpreter','tex');
legend({'Mono','Poly','Bimodal'}, 'Location','best');
box on;

%% ===================== PLOT: l_fc vs beta ==============================
figure('Color','w'); hold on; box on;

for p = 1:nPairs
    if ~isfinite(beta_all(p)) || ~isfinite(lfc_all(p))
        continue;
    end
    
    switch lower(group_all{p})
        case 'mono'
            marker = 's';
            col    = cMono;
        case 'poly'
            marker = 'o';
            col    = cPoly;
        case 'bimo'
            marker = 'd';
            col    = cBimo;
        otherwise
            marker = '^';
            col    = [0 0 0];
    end
    
    plot(beta_all(p), lfc_all(p), marker, ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
end

xlabel('\beta','Interpreter','tex');
ylabel('l_{fc}','Interpreter','tex');   % [nm]
title('Fractocohesive length vs. \beta','Interpreter','tex');
legend({'Mono','Poly','Bimodal'}, 'Location','best');
box on;

fprintf('\nSummary plotting completed.\n');
