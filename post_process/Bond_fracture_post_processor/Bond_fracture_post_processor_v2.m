% ---------------------------------------------------------------
% Post-process LAMMPS bond + atom dumps to compute, for MANY samples:
%   (1) cumulative broken bond count vs macroscopic y-stretch
%   (2) cumulative energy loss due to bond rupture vs macroscopic y-stretch
%   (3) std(N_kuhn) at initial timestep
%
% and then an extra plot:
%   Final accumulated energy loss vs std(N_kuhn) at t = first step
%   - Polydisperse samples: square markers
%   - Monodisperse samples: circle markers
%
% Assumptions:
%   - MATLAB R2016a compatible (no implicit expansion, no "contains", etc.)
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
%   - Bonds only break (no new bonds formed).
%   - Rupture criterion in the simulation: lambda = r / (N*b) > lamc.
%   - Energy per bond at stretch lambda:
%       E_bond(lambda) = N * [ lambda^2/2 - log(1 - lambda^2) ]
%   - We approximate energy loss when a bond breaks as E_bond(lamc),
%     with lamc set to 0.95 (or as desired).
%
%   - We output everything vs macroscopic y-stretch:
%       lambda_y = Ly / Ly0  (Ly0 = box height at first step)
% ---------------------------------------------------------------
clear; clc; close all;

%% ===================== USER SETTINGS =====================

% ---- List of POLYDISPERSE sample directories ----
polyDirs = {
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\sample_alpine_huge_erate_1e-6_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_smp3_notched'
    % add/remove as needed
};

% ---- List of MONODISPERSE sample directories ----
monoDirs = {
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp2_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp3_notched', ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp4_notched' , ...
    'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp5_notched'
    % add/remove as needed
};

% File names used inside each sample directory:
bondDumpName = 'bonds.dump';
atomDumpName = 'atoms1.dump';

% Rupture stretch and clamped region thickness
lamc        = 0.95;   % rupture stretch lambda_c
maskFracY   = 0.12;   % thickness of clamped top/bottom regions (fraction of box height)
debug_print = false;  % per-sample debug printing (set true to inspect)

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
    'finalBroken', {}, ...
    'finalEnergyLoss', {} );

% Pre-allocate vectors for std and final energy
stdN_init_all      = zeros(nSamples,1);
finalEnergy_all    = zeros(nSamples,1);
finalBroken_all    = zeros(nSamples,1);
isPolydisperse_all = false(nSamples,1);

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
    
    % Will store std(N_kuhn) at initial timestep
    stdN_init = NaN;
    
    %% --------- Time-stepping loop for THIS sample ----------
    while true
        % --------- Read one bond timestep ---------
        [okB, tB, Nbonds, xloB, xhiB, yloB, yhiB, zloB, zhiB, ...
            typeB, id1B, id2B, lenB, forceB, NB, bB] = readBondStep(fidB);
        if ~okB
            % End of bond file
            break;
        end
        
        % For this script we only need NB at the first timestep to get std(N_kuhn)
        if stepIndex == 0
            if ~isempty(NB)
                stdN_init = std(double(NB));  % store initial N std
            else
                stdN_init = NaN;
            end
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
        
        stepIndex = stepIndex + 1;
        
        % ---------------- Macroscopic y-stretch ----------------
        Ly = yhiB - yloB;
        if stepIndex == 1
            Ly0 = Ly;   % reference height
            if Ly0 == 0
                error('Initial box height Ly0 is zero for sample %s.', Samples(s).label);
            end
        end
        lambdaY = Ly / Ly0;
        
        % ---------------- Clamp region (this step) -------------
        maskThickness = maskFracY * Ly;
        yClampBottom  = yloB + maskThickness;
        yClampTop     = yhiB - maskThickness;
        
        % ---------------- Atom absolute y positions ------------
        % y = ylo + ys * (yhi - ylo)
        yAbs = yloA + ys .* (yhiA - yloA);
        
        % Build fast map from atom ID -> y coordinate (dense array)
        maxID = max(atomID);
        idxById = zeros(maxID, 1);   % index in atom arrays
        idxById(atomID) = 1:Natoms;
        
        yById = NaN(maxID, 1);
        yById(atomID) = yAbs;
        
        % --------- Precompute interior atom mask (this step) ----
        % interiorMaskThis(id) = true if that atom lies outside grips
        interiorMaskThis = false(maxID, 1);
        validIdx = ~isnan(yById);
        yValid   = yById(validIdx);
        interiorMaskThis(validIdx) = (yValid >= yClampBottom) & (yValid <= yClampTop);
        
        % ------- Compute unique keys for bonds at this step -----
        % Key = type * 1e10 + min(id1,id2)*1e5 + max(id1,id2)
        % to uniquely identify a bond independent of i<->j ordering.
        a1 = min(id1B, id2B);
        a2 = max(id1B, id2B);
        
        % Use int64 to avoid overflow issues
        keysThis = int64(typeB) * int64(1e10) + ...
                   int64(a1)    * int64(1e5)  + ...
                   int64(a2);
        
        % ---------------------------------------------
        % If we have a previous step, detect broken bonds
        % ---------------------------------------------
        if hasPrev
            % Bonds that existed at prev step but not at this step are "lost"
            [~, idxPrevLost] = setdiff(prevKeys, keysThis);
            
            nLostTotal         = numel(idxPrevLost);
            nLostInterior      = 0;
            energyLostThisStep = 0.0;
            
            if nLostTotal > 0
                % Use the precomputed interior mask from previous step
                maskPrev   = prevInteriorMask;
                nMaskPrev  = numel(maskPrev);
                
                for k = 1:nLostTotal
                    idx = idxPrevLost(k);
                    
                    id1_old = prevId1(idx);
                    id2_old = prevId2(idx);
                    N_old   = prevN(idx);
                    
                    % Skip if atom IDs exceed size of previous mask
                    if id1_old > nMaskPrev || id2_old > nMaskPrev
                        continue;
                    end
                    
                    % Check if bond lies in interior at previous step
                    if ~(maskPrev(id1_old) && maskPrev(id2_old))
                        % Ruptured within grips or undefined atoms => ignore
                        continue;
                    end
                    
                    % Energy at rupture: Erupt = N * EruptFactor
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
            % First step: just initialize counts
            if debug_print
                fprintf('Sample %s: t = %d, lambda_y = %.4f: initialization step\n', ...
                    Samples(s).label, tB, lambdaY);
            end
            hasPrev = true;
        end
        
        % Record time series (use bond timestep & macroscopic stretch)
        timeList(stepIndex,1)       = tB;
        lambdaYList(stepIndex,1)    = lambdaY;
        cumBrokenCount(stepIndex,1) = totalBroken;
        cumEnergyLoss(stepIndex,1)  = totalEnergyLoss;
        
        % Update "prev" data for next iteration
        prevKeys         = keysThis;
        prevN            = NB;
        prevId1          = id1B;
        prevId2          = id2B;   % <-- important
        prevYbyId        = yById;
        prevInteriorMask = interiorMaskThis;
        prevYlo          = yloB;
        prevYhi          = yhiB;
    end
    
    % Close files for this sample
    fclose(fidB);
    fclose(fidA);
    
    % Store results for this sample
    Results(s).lambdaY        = lambdaYList;
    Results(s).cumBrokenCount = cumBrokenCount;
    Results(s).cumEnergyLoss  = cumEnergyLoss;
    Results(s).time           = timeList;
    Results(s).stdN_init      = stdN_init;
    Results(s).finalBroken    = totalBroken;
    Results(s).finalEnergyLoss= totalEnergyLoss;
    
    stdN_init_all(s)   = stdN_init;
    finalEnergy_all(s) = totalEnergyLoss;
    finalBroken_all(s) = totalBroken;
    
    fprintf('Sample %s done. Final broken (interior) = %d, final energy loss = %.6g, std(N_init) = %.4g\n', ...
        Samples(s).label, totalBroken, totalEnergyLoss, stdN_init);
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

% --- Cumulative energy loss vs macroscopic y-stretch ---
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
ylabel('Cumulative energy loss (dimensionless, per LAMMPS bond model)');
title(sprintf('Energy loss from bond rupture vs \\lambda_y (\\lambda_c = %.2f)', lamc));
legend({Samples.label}, 'Location', 'best');
box on;

%% ===================== EXTRA PLOT: Energy vs std(N_kuhn) ====

figure; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

% Polydisperse: squares
idxPoly = find(isPolydisperse_all);
if ~isempty(idxPoly)
    plot(stdN_init_all(idxPoly), finalEnergy_all(idxPoly), 's', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
end

% Monodisperse: circles
idxMono = find(~isPolydisperse_all);
if ~isempty(idxMono)
    plot(stdN_init_all(idxMono), finalEnergy_all(idxMono), 'o', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'LineWidth', 1.5);
end

xlabel('Std(N_{kuhn}) at initial timestep');
ylabel('Final accumulated energy loss (interior bonds)');
title('Accumulated rupture energy vs. polydispersity (std N_{kuhn})');

% Build legend entries that say which is which
legendEntries = {};
if ~isempty(idxPoly)
    legendEntries{end+1} = 'Polydisperse (squares)';
end
if ~isempty(idxMono)
    legendEntries{end+1} = 'Monodisperse (circles)';
end
if ~isempty(legendEntries)
    legend(legendEntries, 'Location', 'best');
end
box on;

%% ===================== DONE ================================

fprintf('\nAll samples processed.\n');
for s = 1:nSamples
    fprintf('Sample %-8s: group=%s, std(N_init)=%.4g, final broken=%d, final energy=%.6g\n', ...
        Samples(s).label, Samples(s).group, stdN_init_all(s), finalBroken_all(s), finalEnergy_all(s));
end
