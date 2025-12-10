% ---------------------------------------------------------------
% Multi-sample bond rupture post-processing with Lake–Thomas-style
% normalization and polydispersity measure beta.
%
% For each sample:
%   - Compute cumulative broken bond count vs lambda_y
%   - Compute cumulative energy loss vs lambda_y
%   - Compute std(N_kuhn) at initial timestep
%   - Compute beta from first-timestep N_kuhn distribution:
%         P(N) ~ exp(-N / beta)
%       (monodisperse -> beta = 0)
%   - Compute Lake–Thomas continuum scale T0
%   - Compute fracture energy per crack length:
%         Gdiss = E_diss_total / L_center
%   - Compute normalized dissipation:
%         Enorm = Gdiss / T0
%
% Extra plots:
%   - Enorm vs std(N_kuhn)
%   - E_diss vs std(N_kuhn)
%   - Enorm vs beta
%   - E_diss vs beta
%
% Assumptions:
%   - MATLAB R2016a compatible (no local functions in this script).
%   - Need helper functions:
%       readBondStep.m, readAtomStep.m
%
%   - Bond dump format per timestep:
%       ITEM: TIMESTEP
%       <t>
%       ITEM: NUMBER OF ENTRIES
%       <Nbonds>
%       ITEM: BOX BOUNDS ss ff pp
%       xlo xhi
%       ylo yhi
%       zlo zhi
%       ITEM: ENTRIES c_1[1] c_1[2] c_1[3] c_2[1] c_2[2] c_2[3] c_2[4]
%       type id1 id2 length force N_kuhn b_kuhn
%
%   - Atom dump format per timestep:
%       ITEM: TIMESTEP
%       <t>
%       ITEM: NUMBER OF ATOMS
%       <Natoms>
%       ITEM: BOX BOUNDS ss ff pp
%       xlo xhi
%       ylo yhi
%       zlo zhi
%       ITEM: ATOMS id mol type xs ys zs fx fy fz
%
%   - Bonds only break (no new bonds).
%   - Rupture stretch lambda_c = lamc.
%   - Energy per bond at rupture: N * EruptFactor with
%       EruptFactor = lamc^2/2 - log(1 - lamc^2).
% ---------------------------------------------------------------
clear; clc; close all;

%% ===================== USER SETTINGS =====================

% ---- List of POLYDISPERSE sample directories ----
polyDirs = {
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\sample_alpine_huge_erate_1e-6_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp2_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp3_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp4_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp1_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp2_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp3_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_middist_smp1_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_middist_smp3_notched'
};

% ---- List of MONODISPERSE sample directories ----
monoDirs = {
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp2_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp3_notched', ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp4_notched' , ...
%     'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp5_notched'
};

% File names used inside each sample directory:
bondDumpName = 'bonds.dump';
atomDumpName = 'atoms1.dump';

% Rupture stretch and clamped region thickness
lamc        = 0.95;   % rupture stretch lambda_c
maskFracY   = 0.12;   % thickness of clamped top/bottom regions (fraction of box height)
debug_print = false;  % per-sample debug printing

% --- Centerline & continuum normalization ---
useCenterlineNorm        = false;  % if true, also compute discrete centerline energy
useContinuumNorm         = true;   % use T0-based normalization
LavgMode                 = 'bond'; % 'bond' or 'node' for L_avg
useCenterlineLengthNorm  = true;   % divide E_diss by effective centerline length
centerlineLengthMode     = 'Lx_minus_crack';  % 'Lx_minus_crack' or 'Lx'

% Crack length (same for all samples, in x-units of the box)
crackLength = 1000;  % set to your notch length in same units as x

% --- Beta fit settings ---
nbins_beta = 25;   % number of bins for histogram in beta fit

% ---------------------------------------------------------------
% Precompute energy factor at rupture:
%   E_bond(lambda_c) = N * (lambda_c^2 / 2 - log(1 - lambda_c^2))
% so Erupt_per_bond = N * EruptFactor.
lamc2       = lamc^2;
term1c      = lamc2 / 2.0;
term2c      = log(1.0 - lamc2);
EruptFactor = term1c - term2c;

%% ===================== BUILD SAMPLE LIST =====================

Samples = struct( ...
    'dataDir', {}, ...
    'group'  , {}, ...   % 'poly' or 'mono'
    'label'  , {}, ...
    'bondFile', {}, ...
    'atomFile', {} );

% Polydisperse
for i = 1:numel(polyDirs)
    Samples(end+1).dataDir  = polyDirs{i};
    Samples(end).group      = 'poly';
    Samples(end).label      = sprintf('PD %d', i);
    Samples(end).bondFile   = fullfile(polyDirs{i}, bondDumpName);
    Samples(end).atomFile   = fullfile(polyDirs{i}, atomDumpName);
end

% Monodisperse
for i = 1:numel(monoDirs)
    Samples(end+1).dataDir  = monoDirs{i};
    Samples(end).group      = 'mono';
    Samples(end).label      = sprintf('MD %d', i);
    Samples(end).bondFile   = fullfile(monoDirs{i}, bondDumpName);
    Samples(end).atomFile   = fullfile(monoDirs{i}, atomDumpName);
end

nSamples = numel(Samples);
if nSamples == 0
    error('No samples specified. Please populate polyDirs and/or monoDirs.');
end

%% ===================== STORAGE OVER ALL SAMPLES =============

Results = struct( ...
    'lambdaY', {}, ...
    'cumBrokenCount', {}, ...
    'cumEnergyLoss', {}, ...
    'time', {}, ...
    'stdN_init', {}, ...
    'beta_init', {}, ...
    'finalBroken', {}, ...
    'finalEnergyLoss', {}, ...
    'Ecenter_raw', {}, ...
    'T0', {}, ...
    'L_center', {}, ...
    'Gdiss', {}, ...
    'Enorm', {} );

stdN_init_all      = zeros(nSamples,1);
beta_init_all      = NaN(nSamples,1);
finalEnergy_all    = zeros(nSamples,1);
finalBroken_all    = zeros(nSamples,1);
isPolydisperse_all = false(nSamples,1);
Enorm_all          = NaN(nSamples,1);

%% ===================== MAIN LOOP OVER SAMPLES ===============

for s = 1:nSamples
    
    dataDir      = Samples(s).dataDir;
    bondDumpFile = Samples(s).bondFile;
    atomDumpFile = Samples(s).atomFile;
    isPolydisperse_all(s) = strcmp(Samples(s).group, 'poly');
    
    fprintf('\n============================================\n');
    fprintf('Processing sample %d/%d: %s (%s)\n', ...
        s, nSamples, Samples(s).label, Samples(s).group);
    fprintf('  Bond file: %s\n', bondDumpFile);
    fprintf('  Atom file: %s\n', atomDumpFile);
    
    % ---- Open files ----
    fidB = fopen(bondDumpFile, 'r');
    if fidB == -1
        warning('Could not open bond dump file: %s. Skipping this sample.', bondDumpFile);
        continue;
    end
    
    fidA = fopen(atomDumpFile, 'r');
    if fidA == -1
        warning('Could not open atom dump file: %s. Skipping this sample.', atomDumpFile);
        fclose(fidB);
        continue;
    end
    
    % ---- Per-sample storage/initialization ----
    timeList        = [];
    lambdaYList     = [];
    cumBrokenCount  = [];
    cumEnergyLoss   = [];
    
    totalBroken     = 0;
    totalEnergyLoss = 0;
    
    hasPrev          = false;
    prevKeys         = [];
    prevN            = [];
    prevId1          = [];
    prevId2          = [];
    prevYbyId        = [];
    prevInteriorMask = [];
    prevYlo          = NaN;
    prevYhi          = NaN;
    
    stepIndex        = 0;
    Ly0              = NaN;
    
    T0       = NaN;
    L_center = NaN;
    Gdiss    = NaN;
    Enorm    = NaN;
    
    % std(N_kuhn) and beta at initial timestep
    stdN_init  = NaN;
    beta_local = NaN;
    
    % Centerline energy (discrete, if used)
    energyCenterline = NaN;
    didCenterline    = false;
    
    %% --------- Time-stepping loop for THIS sample ----------
    while true
        % --------- Read one bond timestep ---------
        [okB, tB, Nbonds, xloB, xhiB, yloB, yhiB, zloB, zhiB, ...
            typeB, id1B, id2B, lenB, forceB, NB, bB] = readBondStep(fidB);
        if ~okB
            % End of bond file
            break;
        end
        
        % --------- Read matching atom timestep ---------
        [okA, tA, Natoms, xloA, xhiA, yloA, yhiA, zloA, zhiA, ...
            atomID, xs, ys, zs, fx, fy, fz] = readAtomStep(fidA);
        if ~okA
            warning('Atom dump ended before bond dump for sample %s. Stopping this sample.', Samples(s).label);
            break;
        end
        
        % Sanity check timestep match (not critical)
        if tA ~= tB
            warning('Sample %s: Timestep mismatch: bonds t=%d, atoms t=%d. Using bond timestep.', ...
                Samples(s).label, tB, tA);
        end
        
        % Absolute node coords
        xAbs = xloA + xs .* (xhiA - xloA);
        yAbs = yloA + ys .* (yhiA - yloA);
        zAbs = zloA + zs .* (zhiA - zloA);
        
        stepIndex = stepIndex + 1;
        
        % ---------------- First-timestep statistics ----------------
        if stepIndex == 1
            % ====== A) std(N_kuhn) and beta from N distribution ======
            if ~isempty(NB)
                N_vals = double(NB(:));
                stdN_init = std(N_vals);
                
                N_pos = N_vals(N_vals > 0);
                if ~isempty(N_pos)
                    if std(N_pos) < 1e-8
                        % Monodisperse: all chains same length
                        beta_local = 0;
                    else
                        % Histogram + exponential fit: P(N) ~ exp(-N/beta)
                        Nmin = min(N_pos);
                        Nmax = max(N_pos);
                        if Nmax > Nmin
                            edges   = linspace(Nmin, Nmax, nbins_beta+1);
                            counts  = histcounts(N_pos, edges);
                            centers = 0.5*(edges(1:end-1) + edges(2:end));
                            P       = counts / sum(counts);
                            mask    = (P > 0);
                            xfit    = centers(mask);
                            yfit    = P(mask);
                            if numel(xfit) >= 3
                                logy = log(yfit);
                                pfit = polyfit(xfit, logy, 1);  % log P ~ a + b N
                                bSlope = pfit(1);
                                if bSlope < 0
                                    beta_local = -1 / bSlope;
                                else
                                    beta_local = NaN;
                                end
                            end
                        end
                    end
                else
                    beta_local = NaN;
                end
            else
                stdN_init  = NaN;
                beta_local = NaN;
            end
            
            % ====== B) Continuum T0 and L_center ======
            % B.1 Average N (independent of LavgMode)
            if ~isempty(NB)
                N_avg = mean(double(NB));
            else
                N_avg = NaN;
            end
            
            % B.2 Average length L_avg
            if strcmp(LavgMode, 'bond')
                L_avg = mean(lenB);
            elseif strcmp(LavgMode, 'node')
                % nearest-neighbor distance across nodes
                coords = [xAbs(:), yAbs(:), zAbs(:)];
                Nnodes = size(coords,1);
                minDist = inf(Nnodes,1);
                for in = 1:Nnodes
                    dx = coords(:,1) - coords(in,1);
                    dy = coords(:,2) - coords(in,2);
                    dz = coords(:,3) - coords(in,3);
                    d2 = dx.^2 + dy.^2 + dz.^2;
                    d2(in) = inf;
                    d  = sqrt(d2);
                    minDist(in) = min(d);
                end
                L_avg = mean(minDist);
            else
                error('Unknown LavgMode: %s (use ''bond'' or ''node'')', LavgMode);
            end
            
            % B.3 Geometry & bond density
            Lx0     = xhiB - xloB;
            Ly0_eff = yhiB - yloB;
            area0   = Lx0 * Ly0_eff;
            Nbonds0 = Nbonds;
            if area0 > 0
                rho_b = Nbonds0 / area0;
            else
                rho_b = NaN;
            end
            
            % B.4 Lake–Thomas T0
            if ~isnan(rho_b) && ~isnan(L_avg) && ~isnan(N_avg)
                T0 = 0.5 * rho_b * L_avg * N_avg * EruptFactor;
                if T0 <= 0
                    warning('Sample %s: T0 <= 0, continuum normalization invalid.', Samples(s).label);
                    T0 = NaN;
                end
            else
                warning('Sample %s: Could not compute T0 (NaN in rho_b, L_avg, or N_avg).', Samples(s).label);
                T0 = NaN;
            end
            
            % B.5 Effective centerline length L_center
            if useCenterlineLengthNorm
                if strcmp(centerlineLengthMode, 'Lx_minus_crack')
                    L_center = Lx0 - crackLength;
                else
                    L_center = Lx0;
                end
                if L_center <= 0
                    warning('Sample %s: L_center <= 0, centerline length invalid.', Samples(s).label);
                    L_center = NaN;
                end
            else
                L_center = NaN;
            end
        end
        
        % ---------------- Macroscopic y-stretch ----------------
        Ly = yhiB - yloB;
        if stepIndex == 1
            Ly0 = Ly;
            if Ly0 == 0
                error('Initial box height Ly0 is zero for sample %s.', Samples(s).label);
            end
        end
        lambdaY = Ly / Ly0;
        
        % ---------------- Clamp region (this step) -------------
        maskThickness = maskFracY * Ly;
        yClampBottom  = yloB + maskThickness;
        yClampTop     = yhiB - maskThickness;
        
        % Build fast map from atom ID -> index
        maxID  = max(atomID);
        idxById = zeros(maxID, 1);
        idxById(atomID) = 1:Natoms;
        
        yById = NaN(maxID, 1);
        yById(atomID) = yAbs;
        
        % Interior atom mask
        interiorMaskThis = false(maxID, 1);
        validIdx = ~isnan(yById);
        yValid   = yById(validIdx);
        interiorMaskThis(validIdx) = (yValid >= yClampBottom) & (yValid <= yClampTop);
        
        % ------- Centerline energy at initial timestep ----------
        if useCenterlineNorm && stepIndex == 1 && ~didCenterline
            yMid = 0.5 * (yhiB + yloB);
            energyCenterline = 0.0;
            nCross = 0;
            
            for ib = 1:Nbonds
                id1 = id1B(ib);
                id2 = id2B(ib);
                
                if id1 > maxID || id2 > maxID
                    continue;
                end
                i1 = idxById(id1);
                i2 = idxById(id2);
                if i1 == 0 || i2 == 0
                    continue;
                end
                
                y1 = yAbs(i1);
                y2 = yAbs(i2);
                
                if (y1 - yMid) * (y2 - yMid) <= 0
                    N_local = NB(ib);
                    energyCenterline = energyCenterline + N_local * EruptFactor;
                    nCross = nCross + 1;
                end
            end
            
            if energyCenterline <= 0
                energyCenterline = NaN;
                warning('Sample %s: No valid bonds crossing centerline; centerline energy undefined.', ...
                    Samples(s).label);
            else
                fprintf('Sample %s: centerline @ t0: nCross=%d, E_center=%.6g\n', ...
                    Samples(s).label, nCross, energyCenterline);
            end
            
            didCenterline = true;
        end
        
        % ------- Compute unique keys for bonds at this step ----
        a1 = min(id1B, id2B);
        a2 = max(id1B, id2B);
        
        keysThis = int64(typeB) * int64(1e10) + ...
                   int64(a1)    * int64(1e5)  + ...
                   int64(a2);
        
        % ----------------- Detect broken bonds ------------------
        if hasPrev
            [~, idxPrevLost] = setdiff(prevKeys, keysThis);
            
            nLostTotal         = numel(idxPrevLost);
            nLostInterior      = 0;
            energyLostThisStep = 0.0;
            
            if nLostTotal > 0
                maskPrev   = prevInteriorMask;
                nMaskPrev  = numel(maskPrev);
                
                for kk = 1:nLostTotal
                    idxLost = idxPrevLost(kk);
                    
                    id1_old = prevId1(idxLost);
                    id2_old = prevId2(idxLost);
                    N_old   = prevN(idxLost);
                    
                    if id1_old > nMaskPrev || id2_old > nMaskPrev
                        continue;
                    end
                    
                    if ~(maskPrev(id1_old) && maskPrev(id2_old))
                        continue;
                    end
                    
                    Erupt = N_old * EruptFactor;
                    energyLostThisStep = energyLostThisStep + Erupt;
                    nLostInterior      = nLostInterior + 1;
                end
            end
            
            totalBroken     = totalBroken + nLostInterior;
            totalEnergyLoss = totalEnergyLoss + energyLostThisStep;
            
            if debug_print
                fprintf('Sample %s: t = %d, lambda_y = %.4f: lost = %d, cum = %d, dE = %.3g, E_cum = %.3g\n', ...
                    Samples(s).label, tB, lambdaY, nLostInterior, totalBroken, energyLostThisStep, totalEnergyLoss);
            end
        else
            if debug_print
                fprintf('Sample %s: t = %d, lambda_y = %.4f: initialization step\n', ...
                    Samples(s).label, tB, lambdaY);
            end
            hasPrev = true;
        end
        
        % Record time series
        timeList(stepIndex,1)       = tB;
        lambdaYList(stepIndex,1)    = lambdaY;
        cumBrokenCount(stepIndex,1) = totalBroken;
        cumEnergyLoss(stepIndex,1)  = totalEnergyLoss;
        
        % Update prev data
        prevKeys         = keysThis;
        prevN            = NB;
        prevId1          = id1B;
        prevId2          = id2B;
        prevYbyId        = yById;
        prevInteriorMask = interiorMaskThis;
        prevYlo          = yloB;
        prevYhi          = yhiB;
    end
    
    % Close files for this sample
    fclose(fidB);
    fclose(fidA);
    
    % --- Fracture energy G_diss (energy per unit crack length) ---
    if useCenterlineLengthNorm && ~isnan(L_center) && L_center > 0
        Gdiss = totalEnergyLoss / L_center;
    else
        Gdiss = NaN;
    end
    
    % --- Normalized energy using continuum T0 (energy/area) ---
    if useContinuumNorm && ~isnan(T0) && T0 > 0 && ~isnan(Gdiss)
        Enorm = Gdiss / T0;
    elseif useContinuumNorm && ~isnan(T0) && T0 > 0
        % fallback: normalize by T0 * Lx0 (approx, if needed)
        Enorm = totalEnergyLoss / (T0 * Lx0);
    elseif useCenterlineNorm && ~isnan(energyCenterline) && energyCenterline > 0
        Enorm = totalEnergyLoss / energyCenterline;
    else
        Enorm = NaN;
    end
    
    % Store results for this sample
    Results(s).lambdaY         = lambdaYList;
    Results(s).cumBrokenCount  = cumBrokenCount;
    Results(s).cumEnergyLoss   = cumEnergyLoss;
    Results(s).time            = timeList;
    Results(s).stdN_init       = stdN_init;
    Results(s).beta_init       = beta_local;
    Results(s).finalBroken     = totalBroken;
    Results(s).finalEnergyLoss = totalEnergyLoss;
    Results(s).Ecenter_raw     = energyCenterline;
    Results(s).T0              = T0;
    Results(s).L_center        = L_center;
    Results(s).Gdiss           = Gdiss;
    Results(s).Enorm           = Enorm;
    
    stdN_init_all(s)   = stdN_init;
    beta_init_all(s)   = beta_local;
    finalEnergy_all(s) = totalEnergyLoss;
    finalBroken_all(s) = totalBroken;
    Enorm_all(s)       = Enorm;
    
    fprintf('Sample %s done. Final broken = %d, final energy = %.6g, std(N_init) = %.4g, beta = %.4g, Gdiss = %.4g, Enorm = %.4g\n', ...
        Samples(s).label, totalBroken, totalEnergyLoss, stdN_init, beta_local, Gdiss, Enorm);
end

%% ===================== PLOTS: CURVES =======================

% --- Cumulative broken bond count vs macroscopic y-stretch ---
figure; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
for s = 1:nSamples
    if isempty(Results(s).lambdaY)
        continue;
    end
    if isPolydisperse_all(s)
        ls = '-';   % solid for poly
    else
        ls = '--';  % dashed for mono
    end
    plot(Results(s).lambdaY, Results(s).cumBrokenCount, ls, 'LineWidth', 2);
end
xlabel('\lambda_y (macroscopic stretch in y)');
ylabel('Cumulative broken bond count (interior only)');
title('Broken bonds vs macroscopic y-stretch (all samples)');
legend({Samples.label}, 'Location', 'best');
box on;

% --- Cumulative energy loss vs macroscopic y-stretch (raw) ---
figure; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
for s = 1:nSamples
    if isempty(Results(s).lambdaY)
        continue;
    end
    if isPolydisperse_all(s)
        ls = '-';
    else
        ls = '--';
    end
    plot(Results(s).lambdaY, Results(s).cumEnergyLoss, ls, 'LineWidth', 2);
end
xlabel('\lambda_y (macroscopic stretch in y)');
ylabel('Cumulative energy loss (dimensionless, per bond model)');
title(sprintf('Energy loss from bond rupture vs \\lambda_y (\\lambda_c = %.2f)', lamc));
legend({Samples.label}, 'Location', 'best');
box on;

%% === EXTRA PLOT: Normalized energy vs std(N_kuhn) ===

figure; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

idxPoly = find(isPolydisperse_all);
idxMono = find(~isPolydisperse_all);

% Polydisperse: squares
if ~isempty(idxPoly)
    plot(stdN_init_all(idxPoly), Enorm_all(idxPoly), 's', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
end

% Monodisperse: circles
if ~isempty(idxMono)
    plot(stdN_init_all(idxMono), Enorm_all(idxMono), 'o', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
end

xlabel('Std(N_{kuhn}) at initial timestep');
ylabel('E_{\mathrm{norm}}');
title('Normalized rupture energy vs polydispersity (std N_{kuhn})');

legendEntries = {};
if ~isempty(idxPoly), legendEntries{end+1} = 'Polydisperse (squares)'; end
if ~isempty(idxMono), legendEntries{end+1} = 'Monodisperse (circles)'; end
if ~isempty(legendEntries)
    legend(legendEntries, 'Location', 'best');
end
box on;

%% === EXTRA PLOT: Final energy vs std(N_kuhn) ===

figure; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

if ~isempty(idxPoly)
    plot(stdN_init_all(idxPoly), finalEnergy_all(idxPoly), 's', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
end

if ~isempty(idxMono)
    plot(stdN_init_all(idxMono), finalEnergy_all(idxMono), 'o', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
end

xlabel('Std(N_{kuhn}) at initial timestep');
ylabel('E_{\mathrm{dissip,total}}');
title('Total rupture energy vs polydispersity (std N_{kuhn})');

legendEntries = {};
if ~isempty(idxPoly), legendEntries{end+1} = 'Polydisperse (squares)'; end
if ~isempty(idxMono), legendEntries{end+1} = 'Monodisperse (circles)'; end
if ~isempty(legendEntries)
    legend(legendEntries, 'Location', 'best');
end
box on;

%% === EXTRA PLOT: Final energy vs beta ===

figure; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

if ~isempty(idxPoly)
    maskP = isfinite(beta_init_all(idxPoly)) & isfinite(finalEnergy_all(idxPoly));
    if any(maskP)
        plot(beta_init_all(idxPoly(maskP)), finalEnergy_all(idxPoly(maskP)), 's', ...
            'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
    end
end

if ~isempty(idxMono)
    maskM = isfinite(beta_init_all(idxMono)) & isfinite(finalEnergy_all(idxMono));
    if any(maskM)
        plot(beta_init_all(idxMono(maskM)), finalEnergy_all(idxMono(maskM)), 'o', ...
            'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
    end
end

xlabel('\beta');
ylabel('E_{\mathrm{dissip,total}}');
title('Total rupture energy vs polydispersity \beta');

legendEntries = {};
if ~isempty(idxPoly), legendEntries{end+1} = 'Polydisperse (squares)'; end
if ~isempty(idxMono), legendEntries{end+1} = 'Monodisperse (circles)'; end
if ~isempty(legendEntries)
    legend(legendEntries, 'Location', 'best');
end
box on;

%% === EXTRA PLOT: Normalized energy vs beta ===

figure; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

if ~isempty(idxPoly)
    maskP = isfinite(beta_init_all(idxPoly)) & isfinite(Enorm_all(idxPoly));
    if any(maskP)
        plot(beta_init_all(idxPoly(maskP)), Enorm_all(idxPoly(maskP)), 's', ...
            'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
    end
end

if ~isempty(idxMono)
    maskM = isfinite(beta_init_all(idxMono)) & isfinite(Enorm_all(idxMono));
    if any(maskM)
        plot(beta_init_all(idxMono(maskM)), Enorm_all(idxMono(maskM)), 'o', ...
            'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
    end
end

xlabel('\beta');
ylabel('E_{\mathrm{norm}}');
title('Normalized rupture energy vs polydispersity \beta');

legendEntries = {};
if ~isempty(idxPoly), legendEntries{end+1} = 'Polydisperse (squares)'; end
if ~isempty(idxMono), legendEntries{end+1} = 'Monodisperse (circles)'; end
if ~isempty(legendEntries)
    legend(legendEntries, 'Location', 'best');
end
box on;

%% ===================== DONE (CONSOLE SUMMARY) =================

fprintf('\nAll samples processed.\n');
for s = 1:nSamples
    fprintf('Sample %-8s: group=%s, std(N_init)=%.4g, beta=%.4g, final broken=%d, final energy=%.6g, Enorm=%.4g\n', ...
        Samples(s).label, Samples(s).group, stdN_init_all(s), ...
        beta_init_all(s), finalBroken_all(s), finalEnergy_all(s), Enorm_all(s));
end

%% ===================== SAVE SUMMARIES (LIGHTWEIGHT) =====================

% Global summary table (small)
globalCSV = fullfile(pwd, 'AllSamples_Summary.csv');
fidG = fopen(globalCSV, 'w');
if fidG == -1
    warning('Could not open global summary CSV for writing: %s', globalCSV);
else
    % Header
    fprintf(fidG, 'label,group,stdN_init,beta_init,finalBroken,finalEnergyLoss,T0,L_center,Gdiss,Enorm\n');
    
    for s = 1:nSamples
        
        % Skip empty (e.g., missing files)
        if isempty(Results(s).lambdaY)
            continue;
        end
        
        % ---- Per-sample filenames ----
        outPrefix = fullfile(Samples(s).dataDir, 'rupture_postproc');
        txtFile   = [outPrefix '_summary.txt'];
        csvFile   = [outPrefix '_curves.csv'];
        matFile   = [outPrefix '_summary.mat'];
        
        % ---- 1) Save a small text summary ----
        fidT = fopen(txtFile, 'w');
        if fidT == -1
            warning('Could not open per-sample txt file for writing: %s', txtFile);
        else
            fprintf(fidT, 'Sample label     : %s\n', Samples(s).label);
            fprintf(fidT, 'Group            : %s\n', Samples(s).group);
            fprintf(fidT, 'Data directory   : %s\n', Samples(s).dataDir);
            fprintf(fidT, 'Bond file        : %s\n', Samples(s).bondFile);
            fprintf(fidT, 'Atom file        : %s\n', Samples(s).atomFile);
            fprintf(fidT, '\n');
            fprintf(fidT, 'lamc             : %.6g\n', lamc);
            fprintf(fidT, 'maskFracY        : %.6g\n', maskFracY);
            fprintf(fidT, 'useContinuumNorm : %d\n', useContinuumNorm);
            fprintf(fidT, 'useCenterlineNorm: %d\n', useCenterlineNorm);
            fprintf(fidT, 'useCenterlineLengthNorm: %d\n', useCenterlineLengthNorm);
            fprintf(fidT, 'centerlineLengthMode   : %s\n', centerlineLengthMode);
            fprintf(fidT, 'crackLength      : %.6g\n', crackLength);
            fprintf(fidT, '\n');
            fprintf(fidT, 'std(N_init)      : %.6g\n', Results(s).stdN_init);
            fprintf(fidT, 'beta_init        : %.6g\n', Results(s).beta_init);
            fprintf(fidT, 'finalBroken      : %d\n', Results(s).finalBroken);
            fprintf(fidT, 'finalEnergyLoss  : %.6g\n', Results(s).finalEnergyLoss);
            fprintf(fidT, 'T0               : %.6g\n', Results(s).T0);
            fprintf(fidT, 'L_center         : %.6g\n', Results(s).L_center);
            fprintf(fidT, 'Gdiss            : %.6g\n', Results(s).Gdiss);
            fprintf(fidT, 'Enorm            : %.6g\n', Results(s).Enorm);
            fclose(fidT);
        end
        
        % ---- 2) Save short curves CSV (lambdaY, broken, energy) ----
        fidC = fopen(csvFile, 'w');
        if fidC == -1
            warning('Could not open per-sample csv file for writing: %s', csvFile);
        else
            fprintf(fidC, 'lambdaY,cumBrokenCount,cumEnergyLoss\n');
            nPts = numel(Results(s).lambdaY);
            for k = 1:nPts
                fprintf(fidC, '%.8g,%.8g,%.8g\n', ...
                    Results(s).lambdaY(k), ...
                    Results(s).cumBrokenCount(k), ...
                    Results(s).cumEnergyLoss(k));
            end
            fclose(fidC);
        end
        
        % ---- 3) Save small MAT summary ----
        summaryStruct = struct();
        summaryStruct.label           = Samples(s).label;
        summaryStruct.group           = Samples(s).group;
        summaryStruct.lamc            = lamc;
        summaryStruct.maskFracY       = maskFracY;
        summaryStruct.crackLength     = crackLength;
        summaryStruct.centerlineMode  = centerlineLengthMode;
        summaryStruct.stdN_init       = Results(s).stdN_init;
        summaryStruct.beta_init       = Results(s).beta_init;
        summaryStruct.finalBroken     = Results(s).finalBroken;
        summaryStruct.finalEnergyLoss = Results(s).finalEnergyLoss;
        summaryStruct.T0              = Results(s).T0;
        summaryStruct.L_center        = Results(s).L_center;
        summaryStruct.Gdiss           = Results(s).Gdiss;
        summaryStruct.Enorm           = Results(s).Enorm;
        summaryStruct.lambdaY         = Results(s).lambdaY;
        summaryStruct.cumBrokenCount  = Results(s).cumBrokenCount;
        summaryStruct.cumEnergyLoss   = Results(s).cumEnergyLoss;
        
        try
            save(matFile, '-struct', 'summaryStruct');
        catch ME
            warning('Could not save per-sample MAT file: %s\n%s', matFile, ME.message);
        end
        
        % ---- Add row to GLOBAL summary CSV ----
        fprintf(fidG, '%s,%s,%.8g,%.8g,%d,%.8g,%.8g,%.8g,%.8g,%.8g\n', ...
            Samples(s).label, Samples(s).group, ...
            Results(s).stdN_init, Results(s).beta_init, ...
            Results(s).finalBroken, Results(s).finalEnergyLoss, ...
            Results(s).T0, Results(s).L_center, ...
            Results(s).Gdiss, Results(s).Enorm);
    end
    
    fclose(fidG);
    fprintf('Saved global summary: %s\n', globalCSV);
end

fprintf('Per-sample summaries saved next to each dump directory.\n');
