% ---------------------------------------------------------------
% Post-process LAMMPS bond + atom dumps to compute:
%   (1) cumulative broken bond count vs macroscopic y-stretch
%   (2) cumulative energy loss due to bond rupture vs macroscopic y-stretch
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
%
% ---------------------------------------------------------------
clear; clc; close all;

%% ===================== USER SETTINGS =====================

dataDir = 'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp5_notched';   % change this

bondDumpFile = fullfile(dataDir, 'bonds.dump');
atomDumpFile = fullfile(dataDir, 'atoms1.dump');

lamc          = 0.95;           % rupture stretch lambda_c
maskFracY     = 0.12;           % thickness of clamped top/bottom regions (fraction of box height)
debug_print   = true;          % set to true only when debugging (slows things down)

% ---------------------------------------------------------------
% Precompute energy factor at rupture:
%   E_bond(lambda_c) = N * (lambda_c^2 / 2 - log(1 - lambda_c^2))
% so Erupt_per_bond = N * EruptFactor.
lamc2       = lamc^2;
term1c      = lamc2 / 2.0;
term2c      = log(1.0 - lamc2);
EruptFactor = term1c - term2c;

%% ===================== OPEN FILES =========================

fidB = fopen(bondDumpFile, 'r');
if fidB == -1
    error('Could not open bond dump file: %s', bondDumpFile);
end

fidA = fopen(atomDumpFile, 'r');
if fidA == -1
    fclose(fidB);
    error('Could not open atom dump file: %s', atomDumpFile);
end

%% ===================== STORAGE ============================

timeList           = [];   % timesteps (we keep them just in case)
lambdaYList        = [];   % macroscopic y-stretch Ly/Ly0
cumBrokenCount     = [];   % cumulative # of broken bonds (interior only)
cumEnergyLoss      = [];   % cumulative energy loss (interior only)

totalBroken        = 0;
totalEnergyLoss    = 0;

hasPrev            = false;

prevKeys           = [];
prevN              = [];
prevId1            = [];
prevYbyId          = [];
prevInteriorMask   = [];   % NEW: interior-atom mask from previous step

prevYlo            = NaN;  % kept but no longer needed for mask
prevYhi            = NaN;

stepIndex          = 0;
Ly0                = NaN;  % reference box height (first step)

%% ===================== MAIN LOOP ==========================

while true
    % --------- Read one bond timestep ---------
    [okB, tB, Nbonds, xloB, xhiB, yloB, yhiB, zloB, zhiB, ...
        typeB, id1B, id2B, lenB, forceB, NB, bB] = readBondStep(fidB);
    if ~okB
        % End of file
        break;
    end
    
    % --------- Read matching atom timestep ---------
    [okA, tA, Natoms, xloA, xhiA, yloA, yhiA, zloA, zhiA, ...
        atomID, xs, ys, zs, fx, fy, fz] = readAtomStep(fidA);
    if ~okA
        warning('Atom dump ended before bond dump. Stopping.');
        break;
    end
    
    % Sanity check timestep match
    if tA ~= tB
        warning('Timestep mismatch: bonds t=%d, atoms t=%d. Using bond timestep.', tB, tA);
    end
    
    stepIndex = stepIndex + 1;
    
    % ---------------- Macroscopic y-stretch ----------------
    Ly = yhiB - yloB;
    if stepIndex == 1
        Ly0 = Ly;   % reference height
        if Ly0 == 0
            error('Initial box height Ly0 is zero.');
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
            maskPrev = prevInteriorMask;
            nMaskPrev = numel(maskPrev);
            
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
            fprintf('t = %d, lambda_y = %.4f: lost = %d, cum = %d, dE = %.3g, E_cum = %.3g\n', ...
                tB, lambdaY, nLostInterior, totalBroken, energyLostThisStep, totalEnergyLoss);
        end
    else
        % First step: just initialize counts
        if debug_print
            fprintf('t = %d, lambda_y = %.4f: initialization step\n', tB, lambdaY);
        end
        hasPrev = true;
    end
    
    % Record time series (use bond timestep & macroscopic stretch)
    timeList(stepIndex,1)        = tB;
    lambdaYList(stepIndex,1)     = lambdaY;
    cumBrokenCount(stepIndex,1)  = totalBroken;
    cumEnergyLoss(stepIndex,1)   = totalEnergyLoss;
    
    % Update "prev" data for next iteration
    prevKeys         = keysThis;
    prevN            = NB;
    prevId1          = id1B;
    prevId2          = id2B;   % <<< ADD THIS LINE
    prevYbyId        = yById;
    prevInteriorMask = interiorMaskThis;
    prevYlo          = yloB;
    prevYhi          = yhiB;
end

% Close files
fclose(fidB);
fclose(fidA);

%% ===================== PLOTS ==============================

% Cumulative broken bond count vs macroscopic y-stretch
figure; hold on;
plot(lambdaYList, cumBrokenCount, 'LineWidth', 2);
xlabel('\lambda_y (macroscopic stretch in y)');
ylabel('Cumulative broken bond count (interior only)');
title('Broken bonds vs macroscopic y-stretch');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

% Cumulative energy loss vs macroscopic y-stretch
figure; hold on;
plot(lambdaYList, cumEnergyLoss, 'LineWidth', 2);
xlabel('\lambda_y (macroscopic stretch in y)');
ylabel('Cumulative energy loss (dimensionless, per LAMMPS bond model)');
title(sprintf('Energy loss from bond rupture vs \\lambda_y (\\lambda_c = %.2f)', lamc));
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

% Done.
fprintf('Post-processing complete.\n');
fprintf('Final broken bonds (interior): %d\n', totalBroken);
fprintf('Final cumulative energy loss: %.6g\n', totalEnergyLoss);
