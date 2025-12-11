% -------------------------------------------------------------------------
% Post-processing: paired (unnotched, notched) samples
%   - Supports MONODISPERSE, POLYDISPERSE, BIMODAL groups
%   - For each pair:
%       * Read stress thermo from UNNOTCHED and NOTCHED directories
%       * Plot stress-stretch for unnotched + notched (nice color scheme)
%       * From NOTCHED bond/atom dumps:
%             - Compute beta from N_kuhn distribution at initial step
%             - Compute cumulative bond energy loss E_loss vs lambda_y
%             - Compute Lake-Thomas T0 and Enorm = G_diss / T0
%       * From UNNOTCHED stress:
%             - Compute work of fracture W_frac
%       * Use NOTCHED rupture results:
%             - Compute G_diss = E_loss_total / L_center
%             - Compute fractocohesive length: l_fc = G_diss / W_frac
%
%   - Global plots:
%       * W_frac vs beta
%       * Enorm vs beta
%       * l_fc vs beta
%       (mono / poly / bimodal distinguished by markers/colors)
%
% Requirements:
%   - MATLAB R2016a compatible (script, no local functions).
%   - Helper functions: readBondStep.m, readAtomStep.m
%
% Assumptions:
%   - In each NOTCHED directory:
%       bonds.dump   (bond dump, as in previous scripts)
%       atoms1.dump  (atom dump, as in previous scripts)
%   - In each UNNOTCHED and NOTCHED directory:
%       One thermo file matching pattern 'stress*.txt'
%       with columns:
%         1 Step | 2 v_ntotbonds | 3 c_virial[1] | 4 c_virial[2] | 5 c_virial[4] |
%         6 PotEng | 7 KinEng | 8 v_kper | 9 c_fmax | 10 c_favg | 11 Lx | 12 Ly |
%         13 c_virials[1] | 14 c_virials[2]
% -------------------------------------------------------------------------
clear; clc; close all;

%% ===================== USER INPUT: PAIRS ===============================
% Each row of the following cell arrays is:
%   { 'UNNOTCHED_DIRECTORY_FULL_PATH', 'NOTCHED_DIRECTORY_FULL_PATH' }
%
% Example (NOTE: N x 2 cell array, NOT nested cells):
% monoPairs = {
%     'E:\...\MD_smp1_un', 'E:\...\MD_smp1_notched';
%     'E:\...\MD_smp2_un', 'E:\...\MD_smp2_notched';
% };

% ----------------- MONODISPERSE pairs -----------------
monoPairs = {
   'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp1_longer', ...
   'E:\PhD\My Research\Polydisperse_fracture\PAPER\MD_smp2_notched';
};

% ----------------- POLYDISPERSE pairs -----------------
polyPairs = {
   'E:\PhD\My Research\Polydisperse_fracture\PAPER\sample_alpine_huge_erate_1e-6', ...
   'E:\PhD\My Research\Polydisperse_fracture\PAPER\sample_alpine_huge_erate_1e-6_notched';

   'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp1', ...
   'E:\PhD\My Research\Polydisperse_fracture\PAPER\PD_widerdist_smp1_notched';
};

% ----------------- BIMODAL pairs -----------------
bimoPairs = {
%    'E:\PhD\...\BD_smp1', 'E:\PhD\...\BD_smp1_notched';
};

% ----------------- Common file names / options -----------------
bondDumpName  = 'bonds.dump';
atomDumpName  = 'atoms1.dump';

% Thermo file name pattern INSIDE each directory:
thermoPattern = 'stress*.txt';

% Stress choice: 'virial' -> column 4, 'reduced' -> column 14
stress_choice = 'virial';
% stress_choice = 'reduced';

% Optional downsample for plotting (set =1 to keep all points)
plot_stride = 1;

% Bond breaking / normalization settings (NOTCHED data)
lamc      = 0.95;     % rupture stretch lambda_c
maskFracY = 0.12;     % thickness of clamped regions (fraction of box height)

% Continuum Lake-Thomas normalization
LavgMode                = 'node';            % 'bond' or 'node'
useCenterlineLengthNorm = true;              % divide E_loss by L_center
centerlineLengthMode    = 'Lx_minus_crack';  % 'Lx_minus_crack' or 'Lx'
crackLength             = 1000;              % in same units as x (from box)

% Beta fit settings (from N_kuhn distribution at initial step)
nbins_beta = 25;      % histogram bins for beta fit

% Failure criterion from notched stress:
%   lambda_fail: first stretch where sigma <= frac_drop * sigma_max
frac_drop = 0.4;

%% ===================== PRECOMPUTE ERUPT FACTOR ========================
% E_bond(lambda_c) = N * ( lambda_c^2 / 2 - log(1 - lambda_c^2) )
lamc2       = lamc^2;
term1c      = lamc2 / 2.0;
term2c      = log(1.0 - lamc2);
EruptFactor = term1c - term2c;

%% ===================== BUILD PAIR LIST =================================
Pairs = struct( ...
    'group', {}, ...           % 'mono' | 'poly' | 'bimo'
    'unDir', {}, ...
    'notDir', {}, ...
    'label', {} );

% Monodisperse
for i = 1:size(monoPairs,1)
    Pairs(end+1).group = 'mono';
    Pairs(end).unDir   = monoPairs{i,1};
    Pairs(end).notDir  = monoPairs{i,2};
    Pairs(end).label   = sprintf('Mono_%d', i);
end

% Polydisperse
for i = 1:size(polyPairs,1)
    Pairs(end+1).group = 'poly';
    Pairs(end).unDir   = polyPairs{i,1};
    Pairs(end).notDir  = polyPairs{i,2};
    Pairs(end).label   = sprintf('Poly_%d', i);
end

% Bimodal
for i = 1:size(bimoPairs,1)
    Pairs(end+1).group = 'bimo';
    Pairs(end).unDir   = bimoPairs{i,1};
    Pairs(end).notDir  = bimoPairs{i,2};
    Pairs(end).label   = sprintf('Bimo_%d', i);
end

nPairs = numel(Pairs);
if nPairs == 0
    error('No pairs specified. Please fill monoPairs, polyPairs, and/or bimoPairs.');
end

%% ===================== STORAGE FOR RESULTS =============================
% Store scalars used for summary and global plots
group_all       = cell(nPairs,1);
label_all       = cell(nPairs,1);
beta_all        = nan(nPairs,1);
stdN_all        = nan(nPairs,1);
Wfrac_all       = nan(nPairs,1);
Gdiss_all       = nan(nPairs,1);
Enorm_all       = nan(nPairs,1);
lfail_all       = nan(nPairs,1);
Eloss_total_all = nan(nPairs,1);
T0_all          = nan(nPairs,1);
Lcenter_all     = nan(nPairs,1);
lfc_all         = nan(nPairs,1);   % fractocohesive length

% For stress-stretch plotting
all_lambda_un  = cell(nPairs,1);
all_sigma_un   = cell(nPairs,1);
all_lambda_not = cell(nPairs,1);
all_sigma_not  = cell(nPairs,1);

% For Eloss vs stretch
Eloss_lambda_not = cell(nPairs,1);
Eloss_cum_not    = cell(nPairs,1);

%% ===================== PROCESS EACH PAIR ===============================
fprintf('\n=========== Processing %d pairs ===========\n', nPairs);

for p = 1:nPairs
    group_all{p} = Pairs(p).group;
    label_all{p} = Pairs(p).label;
    
    unDir  = Pairs(p).unDir;
    notDir = Pairs(p).notDir;
    
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Pair %d / %d: %s\n', p, nPairs, Pairs(p).label);
    fprintf('  Unnotched dir: %s\n', unDir);
    fprintf('  Notched   dir: %s\n', notDir);
    
    %% ---- 1) Read thermo (UNNOTCHED + NOTCHED) for stress-stretch ----
    % -------------------- UNNOTCHED --------------------
    d_un = dir(fullfile(unDir, thermoPattern));
    if isempty(d_un)
        error('Pair %d: No thermo file matching pattern "%s" in %s', ...
            p, thermoPattern, unDir);
    end
    thermo_un = fullfile(unDir, d_un(1).name);
    fprintf('  Unnotched thermo: %s\n', thermo_un);
    
    fid = fopen(thermo_un, 'r');
    if fid < 0
        error('Cannot open thermo file: %s', thermo_un);
    end
    
    rows = [];
    tline = fgetl(fid);
    while ischar(tline)
        s = strtrim(tline);
        if ~isempty(s)
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
        tline = fgetl(fid);
    end
    fclose(fid);
    
    if isempty(rows)
        error('No numeric thermo rows found in %s', thermo_un);
    end
    
    Step   = rows(:,1);
    cvir2  = rows(:,4);
    cvirs2 = rows(:,14);
    KinEng = rows(:,7); %#ok<NASGU>
    kper   = rows(:,8); %#ok<NASGU>
    fmax   = rows(:,9); %#ok<NASGU>
    Ly     = rows(:,12);
    
    good = isfinite(Step) & isfinite(Ly) & ...
           (isfinite(cvir2) | isfinite(cvirs2));
    Step   = Step(good);
    cvir2  = cvir2(good);
    cvirs2 = cvirs2(good);
    Ly     = Ly(good);
    
    [~, idxLast] = unique(Step, 'last');
    keep = false(size(Step));
    keep(idxLast) = true;
    
    Step   = Step(keep);
    cvir2  = cvir2(keep);
    cvirs2 = cvirs2(keep);
    Ly     = Ly(keep);
    
    [Step, ord] = sort(Step);
    cvir2  = cvir2(ord);
    cvirs2 = cvirs2(ord);
    Ly     = Ly(ord);
    
    if isempty(Ly) || Ly(1) == 0
        error('Ly(1) is empty or zero for unnotched thermo in pair %d.', p);
    end
    lambda_un = Ly ./ Ly(1);
    
    switch lower(stress_choice)
        case 'virial'
            sigma_un_raw = cvir2;
        case 'reduced'
            sigma_un_raw = cvirs2;
        otherwise
            error('Unknown stress_choice: %s', stress_choice);
    end
    sigma_un = -sigma_un_raw;  % tension positive
    
    if plot_stride > 1
        idx = 1:plot_stride:numel(lambda_un);
        lambda_un = lambda_un(idx);
        sigma_un  = sigma_un(idx);
    end
    
    all_lambda_un{p} = lambda_un;
    all_sigma_un{p}  = sigma_un;
    
    % -------------------- NOTCHED (thermo) --------------------
    d_not = dir(fullfile(notDir, thermoPattern));
    if isempty(d_not)
        error('Pair %d: No thermo file matching pattern "%s" in %s', ...
            p, thermoPattern, notDir);
    end
    thermo_not = fullfile(notDir, d_not(1).name);
    fprintf('  Notched thermo:   %s\n', thermo_not);
    
    fid = fopen(thermo_not, 'r');
    if fid < 0
        error('Cannot open thermo file: %s', thermo_not);
    end
    
    rows = [];
    tline = fgetl(fid);
    while ischar(tline)
        s = strtrim(tline);
        if ~isempty(s)
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
        tline = fgetl(fid);
    end
    fclose(fid);
    
    if isempty(rows)
        error('No numeric thermo rows found in %s', thermo_not);
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
        error('Ly(1) is empty or zero for notched thermo in pair %d.', p);
    end
    lambda_not = Ly_n ./ Ly_n(1);
    
    switch lower(stress_choice)
        case 'virial'
            sigma_not_raw = cvir2_n;
        case 'reduced'
            sigma_not_raw = cvirs2_n;
    end
    sigma_not = -sigma_not_raw;   % tension positive
    
    if plot_stride > 1
        idx = 1:plot_stride:numel(lambda_not);
        lambda_not = lambda_not(idx);
        sigma_not  = sigma_not(idx);
    end
    
    all_lambda_not{p} = lambda_not;
    all_sigma_not{p}  = sigma_not;
    
    %% ---- 2) Failure stretch lambda_fail from NOTCHED stress ----
    sig = sigma_not;
    lam = lambda_not;
    nloc = min(numel(sig), numel(lam));
    sig = sig(1:nloc);
    lam = lam(1:nloc);
    
    [sig_max, imax] = max(sig);
    lambda_fail = NaN;
    if ~isempty(imax) && isfinite(sig_max)
        sig_post = sig(imax:end);
        lam_post = lam(imax:end);
        thresh   = frac_drop * sig_max;
        idx_drop = find(sig_post <= thresh, 1, 'first');
        if ~isempty(idx_drop)
            lambda_fail = lam_post(idx_drop);
        else
            lambda_fail = lam(imax);
        end
    end
    lfail_all(p) = lambda_fail;
    fprintf('  lambda_fail (notched) = %.4f\n', lambda_fail);
    
    %% ---- 3) Work of fracture W_frac from UNNOTCHED stress ----
    Wfrac = NaN;
    if ~isnan(lambda_fail)
        lu  = lambda_un;
        su  = sigma_un;
        mask = (lu >= 1) & (lu <= lambda_fail);
        lam_clip = lu(mask);
        sig_clip = su(mask);
        if numel(lam_clip) >= 2
            Wfrac = trapz(lam_clip, sig_clip);
        else
            warning('Pair %d: Not enough points between 1 and lambda_fail for W_frac.', p);
        end
    else
        warning('Pair %d: lambda_fail is NaN, W_frac undefined.', p);
    end
    Wfrac_all(p) = Wfrac;
    fprintf('  W_frac = %.4e\n', Wfrac);
    
    %% ---- 4) Rupture post-processing from NOTCHED bond/atom dumps ----
    bondDumpFile = fullfile(notDir, bondDumpName);
    atomDumpFile = fullfile(notDir, atomDumpName);
    
    fprintf('  Notched bond dump: %s\n', bondDumpFile);
    fprintf('  Notched atom dump: %s\n', atomDumpFile);
    
    fidB = fopen(bondDumpFile, 'r');
    if fidB == -1
        warning('Pair %d: Could not open bond dump; skipping rupture-based quantities.', p);
        continue;
    end
    
    fidA = fopen(atomDumpFile, 'r');
    if fidA == -1
        warning('Pair %d: Could not open atom dump; skipping rupture-based quantities.', p);
        fclose(fidB);
        continue;
    end
    
    % Per-pair storage
    timeList       = [];
    lambdaYList    = [];
    cumEnergyLoss  = [];
    
    totalEnergyLoss = 0;
    stepIndex       = 0;
    Ly0_bond        = NaN;
    
    hasPrev          = false;
    prevKeys         = [];
    prevN            = [];
    prevId1          = [];
    prevId2          = [];
    prevInteriorMask = [];
    
    stdN_init  = NaN;
    beta_local = NaN;
    T0         = NaN;
    L_center   = NaN;
    
    % Time-stepping over bond/atom dumps
    while true
        [okB, tB, Nbonds, xloB, xhiB, yloB, yhiB, zloB, zhiB, ...
            typeB, id1B, id2B, lenB, forceB, NB, bB] = readBondStep(fidB); %#ok<NASGU>
        if ~okB
            break;
        end
        
        [okA, tA, Natoms, xloA, xhiA, yloA, yhiA, zloA, zhiA, ...
            atomID, xs, ys, zs, fxA, fyA, fzA] = readAtomStep(fidA); %#ok<NASGU>
        if ~okA
            warning('Pair %d: Atom dump ended before bond dump.', p);
            break;
        end
        
        xAbs = xloA + xs .* (xhiA - xloA);
        yAbs = yloA + ys .* (yhiA - yloA);
        zAbs = zloA + zs .* (zhiA - zloA); %#ok<NASGU>
        
        stepIndex = stepIndex + 1;
        
        % ----- Initial timestep: beta, std(N), T0, L_center -----
        if stepIndex == 1
            if ~isempty(NB)
                N_vals = double(NB(:));
                stdN_init = std(N_vals);
                
                N_pos = N_vals(N_vals > 0);
                if ~isempty(N_pos)
                    if std(N_pos) < 1e-8
                        % effectively monodisperse
                        if strcmp(Pairs(p).group, 'mono')
                            beta_local = 0;
                        else
                            beta_local = 0; % or near-0
                        end
                    else
                        Nmin = min(N_pos);
                        Nmax = max(N_pos);
                        if Nmax > Nmin
                            edges   = linspace(Nmin, Nmax, nbins_beta+1);
                            counts  = histcounts(N_pos, edges);
                            centers = 0.5*(edges(1:end-1) + edges(2:end));
                            P       = counts / sum(counts);
                            maskP   = (P > 0);
                            xfit    = centers(maskP);
                            yfit    = P(maskP);
                            if numel(xfit) >= 3
                                logy = log(yfit);
                                pf   = polyfit(xfit, logy, 1);   % log P ~ a + b N
                                b_slope = pf(1);
                                if b_slope < 0
                                    beta_local = -1 / b_slope;
                                else
                                    beta_local = NaN;
                                end
                            end
                        end
                    end
                end
            end
            
            % Override for monodisperse: by definition beta = 0
            if strcmp(Pairs(p).group, 'mono')
                beta_local = 0;
            end
            
            % --- Lake-Thomas T0 and L_center ---
            if ~isempty(NB)
                N_avg = mean(double(NB));
            else
                N_avg = NaN;
            end
            
            if strcmp(LavgMode, 'bond')
                L_avg = mean(lenB);
            elseif strcmp(LavgMode, 'node')
                coords = [xAbs(:), yAbs(:), zAbs(:)];
                Nnodes = size(coords,1);
                minDist = inf(Nnodes,1);
                for in = 1:Nnodes
                    dx = coords(:,1) - coords(in,1);
                    dy = coords(:,2) - coords(in,2);
                    dz = coords(:,3) - coords(in,3);
                    d2 = dx.^2 + dy.^2 + dz.^2;
                    d2(in) = inf;
                    d = sqrt(d2);
                    minDist(in) = min(d);
                end
                L_avg = mean(minDist);
            else
                error('Unknown LavgMode: %s', LavgMode);
            end
            
            Lx0     = xhiB - xloB;
            Ly0_eff = yhiB - yloB;
            area0   = Lx0 * Ly0_eff;
            Nbonds0 = Nbonds;
            if area0 > 0
                rho_b = Nbonds0 / area0;
            else
                rho_b = NaN;
            end
            
            if ~isnan(rho_b) && ~isnan(L_avg) && ~isnan(N_avg)
                T0 = 0.5 * rho_b * L_avg * N_avg * EruptFactor;
                if T0 <= 0
                    warning('Pair %d: T0 <= 0, continuum normalization invalid.', p);
                    T0 = NaN;
                end
            else
                warning('Pair %d: Could not compute T0 (NaN in rho_b, L_avg, or N_avg).', p);
                T0 = NaN;
            end
            
            if useCenterlineLengthNorm
                if strcmp(centerlineLengthMode, 'Lx_minus_crack')
                    L_center = Lx0 - crackLength;
                else
                    L_center = Lx0;
                end
                if L_center <= 0
                    warning('Pair %d: L_center <= 0.', p);
                    L_center = NaN;
                end
            else
                L_center = NaN;
            end
        end
        
        % --- Macroscopic stretch from current bond box ---
        Ly_bond = yhiB - yloB;
        if stepIndex == 1
            Ly0_bond = Ly_bond;
            if Ly0_bond == 0
                error('Pair %d: Ly0_bond is zero.', p);
            end
        end
        lambdaY = Ly_bond / Ly0_bond;
        
        % --- Interior mask this step ---
        maskThickness = maskFracY * Ly_bond;
        yClampBottom  = yloB + maskThickness;
        yClampTop     = yhiB - maskThickness;
        
        maxID  = max(atomID);
        idxById = zeros(maxID,1);
        idxById(atomID) = 1:Natoms;
        
        yById = NaN(maxID,1);
        yById(atomID) = yAbs;
        
        interiorMaskThis = false(maxID,1);
        validIdx = ~isnan(yById);
        yValid   = yById(validIdx);
        interiorMaskThis(validIdx) = (yValid >= yClampBottom) & (yValid <= yClampTop);
        
        % --- Bond keys for this step ---
        a1 = min(id1B, id2B);
        a2 = max(id1B, id2B);
        keysThis = int64(typeB) * int64(1e10) + ...
                   int64(a1)    * int64(1e5)  + ...
                   int64(a2);
        
        % --- Bond energy loss from broken interior bonds ---
        if hasPrev
            [~, idxPrevLost] = setdiff(prevKeys, keysThis);
            energyLostThisStep = 0.0;
            nMaskPrev = numel(prevInteriorMask);
            
            for kk = 1:numel(idxPrevLost)
                idxLost = idxPrevLost(kk);
                id1_old = prevId1(idxLost);
                id2_old = prevId2(idxLost);
                N_old   = prevN(idxLost);
                
                if id1_old > nMaskPrev || id2_old > nMaskPrev
                    continue;
                end
                if ~(prevInteriorMask(id1_old) && prevInteriorMask(id2_old))
                    continue;
                end
                
                Erupt = N_old * EruptFactor;
                energyLostThisStep = energyLostThisStep + Erupt;
            end
            
            totalEnergyLoss = totalEnergyLoss + energyLostThisStep;
        else
            hasPrev = true;
        end
        
        timeList(stepIndex,1)      = tB;
        lambdaYList(stepIndex,1)   = lambdaY;
        cumEnergyLoss(stepIndex,1) = totalEnergyLoss;
        
        prevKeys         = keysThis;
        prevN            = NB;
        prevId1          = id1B;
        prevId2          = id2B;
        prevInteriorMask = interiorMaskThis;
    end
    
    fclose(fidB);
    fclose(fidA);
    
    Eloss_total = totalEnergyLoss;
    Eloss_total_all(p) = Eloss_total;
    
    % G_diss = totalEnergyLoss / L_center
    if useCenterlineLengthNorm && ~isnan(L_center) && L_center > 0
        Gdiss = Eloss_total / L_center;
    else
        Gdiss = NaN;
    end
    Gdiss_all(p) = Gdiss;
    
    % Enorm = G_diss / T0
    Enorm = NaN;
    if ~isnan(Gdiss) && ~isnan(T0) && T0 > 0
        Enorm = Gdiss / T0;
    end
    Enorm_all(p) = Enorm;
    
    % Fractocohesive length: l_fc = G_diss / W_frac
    lfc = NaN;
    if ~isnan(Gdiss) && Gdiss > 0 && ~isnan(Wfrac) && Wfrac > 0
        lfc = Gdiss / Wfrac;
    end
    lfc_all(p) = lfc;
    
    % Store local scalars / curves
    beta_all(p)         = beta_local;
    stdN_all(p)         = stdN_init;
    T0_all(p)           = T0;
    Lcenter_all(p)      = L_center;
    Eloss_lambda_not{p} = lambdaYList;
    Eloss_cum_not{p}    = cumEnergyLoss;
    
    fprintf('  std(N_init) = %.4g, beta = %.4g\n', stdN_init, beta_local);
    fprintf('  Eloss_total = %.6g, Gdiss = %.6g, Enorm = %.6g, l_fc = %.6g\n', ...
        Eloss_total, Gdiss, Enorm, lfc);
    
    %% ===================== NEW: PER-PAIR TXT SUMMARY ====================
    pairSummaryTxt = fullfile(notDir, [Pairs(p).label '_pair_summary.txt']);
    fidP = fopen(pairSummaryTxt, 'w');
    if fidP == -1
        warning('Pair %d: Could not open per-pair summary txt for writing: %s', ...
            p, pairSummaryTxt);
    else
        fprintf(fidP, 'Pair label       : %s\n', Pairs(p).label);
        fprintf(fidP, 'Group            : %s\n', Pairs(p).group);
        fprintf(fidP, 'Unnotched dir    : %s\n', unDir);
        fprintf(fidP, 'Notched dir      : %s\n', notDir);

        fprintf(fidP, '\n--- Global settings ---\n');
        fprintf(fidP, 'lamc             : %.6g      [-]\n', lamc);
        fprintf(fidP, 'maskFracY        : %.6g      [-]\n', maskFracY);
        fprintf(fidP, 'LavgMode         : %s\n', LavgMode);
        fprintf(fidP, 'useCenterlineLengthNorm : %d\n', useCenterlineLengthNorm);
        fprintf(fidP, 'centerlineLengthMode    : %s\n', centerlineLengthMode);
        fprintf(fidP, 'crackLength      : %.6g      [nm]\n', crackLength);

        fprintf(fidP, '\n--- Rupture / polydispersity ---\n');
        fprintf(fidP, 'std(N_init)      : %.6g      [-]\n', stdN_init);
        fprintf(fidP, 'beta             : %.6g      [-]\n', beta_local);

        fprintf(fidP, '\n--- Energetic quantities ---\n');
        fprintf(fidP, 'W_frac           : %.6g      [MPa]\n', Wfrac);
        fprintf(fidP, 'Eloss_total      : %.6g      [1e-21 J]\n', Eloss_total);
        fprintf(fidP, 'T0               : %.6g      [1e-21 J/nm^2]\n', T0);
        fprintf(fidP, 'L_center         : %.6g      [nm]\n', L_center);
        fprintf(fidP, 'Gdiss            : %.6g      [1e-21 J/nm^2]\n', Gdiss);
        fprintf(fidP, 'Enorm            : %.6g      [-]\n', Enorm);
        fprintf(fidP, 'lambda_fail      : %.6g      [-]\n', lambda_fail);
        fprintf(fidP, 'l_fc             : %.6g      [nm]\n', lfc);

        fclose(fidP);
        fprintf('  Saved per-pair summary txt: %s\n', pairSummaryTxt);
    end

    %% ====================================================================
end

%% ===================== PLOTTING: STRESS–STRETCH =======================
figsize = [1 1 4 3];
set(0,'DefaultFigureUnits','inches','DefaultFigurePosition',figsize);
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',14);

cMono = [203 54 88]/255;
cPoly = [22 114 121]/255;
cBimo = [0.2 0.2 0.2];

figure('Color','w'); hold on; box on;

for p = 1:nPairs
    lam_un  = all_lambda_un{p};
    sig_un  = all_sigma_un{p};
    lam_not = all_lambda_not{p};
    sig_not = all_sigma_not{p};
    
    if isempty(lam_un) || isempty(sig_un) || isempty(lam_not) || isempty(sig_not)
        continue;
    end
    
    switch group_all{p}
        case 'mono'
            col = cMono;
        case 'poly'
            col = cPoly;
        case 'bimo'
            col = cBimo;
        otherwise
            col = [0 0 0];
    end
    
    plot(lam_un,  sig_un,  '-',  'LineWidth', 1.5, 'Color', col); % unnotched
    plot(lam_not, sig_not, '--', 'LineWidth', 1.5, 'Color', col); % notched
end

xlabel('\lambda_y','Interpreter','tex');
ylabel('\sigma_{yy}','Interpreter','tex');
title('\sigma_{yy} vs. \lambda_y (unnotched solid, notched dashed)');
box on;

%% ===================== PLOTTING: W_frac vs beta =======================
figure('Color','w'); hold on; box on;

for p = 1:nPairs
    if ~isfinite(beta_all(p)) || ~isfinite(Wfrac_all(p))
        continue;
    end
    switch group_all{p}
        case 'mono'
            marker = 's';
            col = cMono;
        case 'poly'
            marker = 'o';
            col = cPoly;
        case 'bimo'
            marker = 'd';
            col = cBimo;
        otherwise
            marker = '^';
            col = [0 0 0];
    end
    plot(beta_all(p), Wfrac_all(p), marker, ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
end

xlabel('\beta','Interpreter','tex');
ylabel('W_{frac}','Interpreter','tex');
title('Work of fracture vs. \beta','Interpreter','tex');
legend({'Mono','Poly','Bimodal'}, 'Location','best');
box on;



%% ===================== PLOTTING: Enorm vs beta ========================
figure('Color','w'); hold on; box on;

for p = 1:nPairs
    if ~isfinite(beta_all(p)) || ~isfinite(Enorm_all(p))
        continue;
    end
    switch group_all{p}
        case 'mono'
            marker = 's';
            col = cMono;
        case 'poly'
            marker = 'o';
            col = cPoly;
        case 'bimo'
            marker = 'd';
            col = cBimo;
        otherwise
            marker = '^';
            col = [0 0 0];
    end
    plot(beta_all(p), Enorm_all(p), marker, ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
end

xlabel('\beta','Interpreter','tex');
ylabel('E_{norm} = G_{diss}/T_0','Interpreter','tex');
title('Normalized rupture energy vs. \beta','Interpreter','tex');
legend({'Mono','Poly','Bimodal'}, 'Location','best');
box on;


%% ===================== PLOTTING: l_fc vs beta =========================
figure('Color','w'); hold on; box on;

for p = 1:nPairs
    if ~isfinite(beta_all(p)) || ~isfinite(lfc_all(p))
        continue;
    end
    switch group_all{p}
        case 'mono'
            marker = 's';
            col = cMono;
        case 'poly'
            marker = 'o';
            col = cPoly;
        case 'bimo'
            marker = 'd';
            col = cBimo;
        otherwise
            marker = '^';
            col = [0 0 0];
    end
    plot(beta_all(p), lfc_all(p), marker, ...
        'MarkerSize', 8, 'LineWidth', 1.5, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', col);
end

xlabel('\beta','Interpreter','tex');
ylabel('l_{fc} = G_{diss}/W_{frac}','Interpreter','tex');
title('Fractocohesive length vs. \beta','Interpreter','tex');
legend({'Mono','Poly','Bimodal'}, 'Location','best');
box on;


%% ===================== SAVE GLOBAL SUMMARY ============================
summaryCSV = fullfile(pwd, 'Pairs_Summary.csv');
fidG = fopen(summaryCSV, 'w');
if fidG == -1
    warning('Could not open global summary CSV for writing: %s', summaryCSV);
else
    % Header
    fprintf(fidG, ['pair_label,group,beta,stdN_init,Wfrac,Eloss_total,' ...
                   'Gdiss,Enorm,lfc,lambda_fail,T0,L_center\n']);
    
    for p = 1:nPairs
        if isnan(beta_all(p)) && isnan(Wfrac_all(p)) && isnan(Enorm_all(p))
            continue;
        end
        fprintf(fidG, '%s,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n', ...
            label_all{p}, group_all{p}, ...
            beta_all(p), stdN_all(p), Wfrac_all(p), Eloss_total_all(p), ...
            Gdiss_all(p), Enorm_all(p), lfc_all(p), lfail_all(p), ...
            T0_all(p), Lcenter_all(p));
    end
    fclose(fidG);
    fprintf('\nSaved global summary CSV: %s\n', summaryCSV);
end

% Also save a MAT with the same scalars + some curves
PairsSummary = struct();
PairsSummary.group_all        = group_all;
PairsSummary.label_all        = label_all;
PairsSummary.beta_all         = beta_all;
PairsSummary.stdN_all         = stdN_all;
PairsSummary.Wfrac_all        = Wfrac_all;
PairsSummary.Eloss_total_all  = Eloss_total_all;
PairsSummary.Gdiss_all        = Gdiss_all;
PairsSummary.Enorm_all        = Enorm_all;
PairsSummary.lfc_all          = lfc_all;
PairsSummary.lfail_all        = lfail_all;
PairsSummary.T0_all           = T0_all;
PairsSummary.Lcenter_all      = Lcenter_all;
PairsSummary.all_lambda_un    = all_lambda_un;
PairsSummary.all_sigma_un     = all_sigma_un;
PairsSummary.all_lambda_not   = all_lambda_not;
PairsSummary.all_sigma_not    = all_sigma_not;
PairsSummary.Eloss_lambda_not = Eloss_lambda_not;
PairsSummary.Eloss_cum_not    = Eloss_cum_not;

summaryMAT = fullfile(pwd, 'Pairs_Summary.mat');
try
    save(summaryMAT, '-struct', 'PairsSummary');
    fprintf('Saved global summary MAT: %s\n', summaryMAT);
catch ME
    warning('Could not save global summary MAT: %s\n%s', summaryMAT, ME.message);
end

fprintf('\nAll pairs processed. Summary files written in current folder.\n');
