% ============================================================
% Plot_RuptureSummaries.m
% ------------------------------------------------------------
% Loads the *_summary.mat and *_curves.csv files produced by
% the main rupture post-processing script.
%
% Then plots:
%   - Broken bonds vs lambda_y
%   - Energy loss vs lambda_y
%   - Normalized energy vs std(N_kuhn)
%
% No LAMMPS dumps are read here. Very fast.
%
% MATLAB R2016a compatible.
% ============================================================

clear; clc; close all;

%% ============== USER SETTINGS ==============================
% Root folders that contain sample directories
rootDirs = {
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp3_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp4_notched' , ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp5_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\sample_alpine_huge_erate_1e-6_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp3_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp4_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp1_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp3_notched',...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_middist_smp1_notched',...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_middist_smp2_notched',...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_middist_smp3_notched'
};
% Add/remove folders as needed.

% Filenames created earlier:
summaryMatName = 'rupture_postproc_summary.mat';
curvesCSVName  = 'rupture_postproc_curves.csv';

% ---------- COLOR / STYLE KNOBS (edit here) -----------------
% Colors for mono (MD) and poly samples:
colorMono = [203 54 88]/255;       % black
colorPoly = [22  114 121]/255;    % dark red

% Line styles used to encode groups of similar std(N_kuhn):
stdTol   = 1.0;            % samples within +/- stdTol share linestyle
lsList   = {'-','--',':','-.'};  % cycles if more groups than entries
% ------------------------------------------------------------

%% ============================================================

Samples = struct('folder', {}, 'label', {}, 'group', {}, ...
                 'lambdaY', {}, 'broken', {}, 'energy', {}, ...
                 'stdN', {}, 'Enorm', {});

for r = 1:numel(rootDirs)
    folder = rootDirs{r};
    
    matFile = fullfile(folder, summaryMatName);
    csvFile = fullfile(folder, curvesCSVName);

    % Skip missing
    if ~exist(matFile, 'file') || ~exist(csvFile, 'file')
        fprintf('Skipping folder (no summary files): %s\n', folder);
        continue;
    end
    
    % -------- Load summary MAT file (small) --------
    S = load(matFile);
    
    % -------- Load curves CSV --------
    fid = fopen(csvFile, 'r');
    if fid == -1
        warning('Cannot open CSV: %s', csvFile);
        continue;
    end
    
    % Skip header
    header = fgetl(fid); %#ok<NASGU>
    
    lambdaY = [];
    broken  = [];
    energyL = [];
    
    line = fgetl(fid);
    while ischar(line)
        values = sscanf(line, '%f,%f,%f');
        if numel(values) == 3
            lambdaY(end+1,1) = values(1);
            broken(end+1,1)  = values(2);
            energyL(end+1,1) = values(3);
        end
        line = fgetl(fid);
    end
    fclose(fid);
    
    % Store result
    idx = numel(Samples)+1;
    Samples(idx).folder = folder;
    Samples(idx).label  = S.label;
    Samples(idx).group  = S.group;   % 'mono' or 'poly'
    Samples(idx).lambdaY = lambdaY;
    Samples(idx).broken  = broken;
    Samples(idx).energy  = energyL;
    Samples(idx).stdN    = S.stdN_init;
    Samples(idx).Enorm   = S.Enorm;
end

nSamples = numel(Samples);
if nSamples == 0
    error('No valid summary files found.');
end

%% ======== PRECOMPUTE STD-BASED LINESTYLE GROUPS =============

stdVals = zeros(nSamples,1);
for s = 1:nSamples
    stdVals(s) = Samples(s).stdN;
end

clusterCenters = [];
clusterIdx = zeros(nSamples,1);

for s = 1:nSamples
    val = stdVals(s);
    assigned = false;
    for c = 1:numel(clusterCenters)
        if abs(val - clusterCenters(c)) <= stdTol
            clusterIdx(s) = c;
            assigned = true;
            break;
        end
    end
    if ~assigned
        clusterCenters(end+1) = val; %#ok<AGROW>
        clusterIdx(s) = numel(clusterCenters);
    end
end

% Map each cluster to a linestyle from lsList
clusterLS = cell(numel(clusterCenters),1);
for c = 1:numel(clusterCenters)
    idxLS = mod(c-1, numel(lsList)) + 1;
    clusterLS{c} = lsList{idxLS};
end

%% ======================= PLOTS ==============================

% -------- Broken bonds vs lambda_y --------
figure; hold on;
set(gca,'FontName','Times New Roman','FontSize',14);

for s = 1:nSamples
    if strcmp(Samples(s).group,'poly')
        col = colorPoly;
        lsGroup = '-';
    else
        col = colorMono;
        lsGroup = '--';
    end
    
    plot(Samples(s).lambdaY, Samples(s).broken, ...
         'LineWidth', 2, 'Color', col, 'LineStyle', lsGroup);
end

xlabel('\lambda_y');
ylabel('Cumulative broken bonds');
title('Broken bonds vs stretch (loaded from summary files)');
legend({Samples.label}, 'Location', 'best');
box on;

% -------- Energy loss vs lambda_y --------
figure; hold on;
set(gca,'FontName','Times New Roman','FontSize',14);

legStr = cell(nSamples,1);

for s = 1:nSamples
    % Color by group
    if strcmp(Samples(s).group,'poly')
        col = colorPoly;
    else
        col = colorMono;
    end
    
    % Line style by std(N_kuhn) cluster
    cIdx = clusterIdx(s);
    lsStd = clusterLS{cIdx};
    
    plot(Samples(s).lambdaY, Samples(s).energy, ...
         'LineWidth', 2, 'Color', col, 'LineStyle', lsStd);
     
    legStr{s} = sprintf('%s  (stdN = %.2f)', ...
                        Samples(s).label, Samples(s).stdN);
end

xlabel('\lambda_y');
ylabel('Energy loss (cumulative)');
title('Energy loss vs stretch (std(N_{kuhn}) encoded in linestyle)');
legend(legStr, 'Location', 'best');
box on;

% -------- Normalized energy vs std(N_kuhn) --------
figure; hold on;
set(gca,'FontName','Times New Roman','FontSize',14);

% scatter filled for poly samples (squares)
for s = 1:nSamples
    if strcmp(Samples(s).group,'poly')
        h = scatter(Samples(s).stdN, Samples(s).Enorm, ...
                    70, 'MarkerFaceColor', colorPoly, ...
                    'MarkerEdgeColor', 'k', ...
                    'Marker', 's', ...
                    'LineWidth', 1.2);
    end
end

% scatter filled for mono samples (circles)
for s = 1:nSamples
    if strcmp(Samples(s).group,'mono')
        h = scatter(Samples(s).stdN, Samples(s).Enorm, ...
                    70, 'MarkerFaceColor', colorMono, ...
                    'MarkerEdgeColor', 'k', ...
                    'Marker', 'o', ...
                    'LineWidth', 1.2);
    end
end

% Dummy handles for legend
hPoly = scatter(NaN,NaN,70, 'MarkerFaceColor', colorPoly, ...
                'MarkerEdgeColor','k','Marker','s','LineWidth',1.2);
hMono = scatter(NaN,NaN,70, 'MarkerFaceColor', colorMono, ...
                'MarkerEdgeColor','k','Marker','o','LineWidth',1.2);

% Dummy handles for clean legend
hPoly = scatter(NaN,NaN,70,colorPoly,'filled'); set(hPoly,'Marker','s');
hMono = scatter(NaN,NaN,70,colorMono,'filled'); set(hMono,'Marker','o');

xlabel('std(N_{kuhn})');
ylabel('Normalized dissipated energy');
title('Enorm vs polydispersity (from summary files)');
legend([hPoly hMono], {'Polydisperse','Monodisperse'}, 'Location','best');
box on;

fprintf('\nLoaded and plotted %d samples.\n', nSamples);
