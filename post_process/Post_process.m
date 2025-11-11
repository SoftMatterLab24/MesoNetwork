% PlotBondDistributions.m
% MATLAB R2016a compatible
% Reads:
%   - bonds.dump  (columns: bond_type, batom1, batom2, end_to_end_length, force)
%   - bond.table  (columns: index,   batom1, batom2, N_kuhn, kuhn_length_b)
%
% Produces:
%   1) Histogram: end-to-end length r
%   2) Histogram: contour length in Kuhn segments N
%   3) Histogram: pre-stretch lambda0 = r/(N*b)

clear; clc;
close all
%% ----------- User files -----------
dataDir = 'E:\PhD\My Research\Polydisperse_fracture\PAPER\sample_alpine_huge_erate_1e-6';
bondsDumpFile = fullfile(dataDir,'bonds.dump');
bondTableFile = fullfile(dataDir,'bond.table');

assert(exist(bondsDumpFile,'file')==2, 'Missing: %s', bondsDumpFile);
assert(exist(bondTableFile,'file')==2, 'Missing: %s', bondTableFile);
fprintf('Using:\n  %s\n  %s\n', bondsDumpFile, bondTableFile);

%% ----------- Plot mode knob -----------
% 'hist' (default, your current look) or 'line' (P(x) vs x)
plotMode = 'hist';  % 'hist' or 'line'


%% ----------- Parse bonds.dump (only timestep 0) -----------
fprintf('Reading first (t=0) block from %s ...\n', bondsDumpFile);
fid = fopen(bondsDumpFile,'r');
assert(fid>0,'Cannot open %s', bondsDumpFile);

B = []; reading = false; done = false;
while ~feof(fid) && ~done
    ln = fgetl(fid);
    if isempty(ln), continue; end

    if strncmp(ln,'ITEM: TIMESTEP',14)
        % next line is timestep value
        ts = str2double(fgetl(fid));
        if ts == 0
            reading = true;
        else
            if reading
                done = true; % we already passed t=0 block
            end
        end
        continue;
    end

    if reading && strncmp(ln,'ITEM: ENTRIES',13)
        % start reading numeric block
        while true
            pos = ftell(fid);
            L = fgetl(fid);
            if ~ischar(L) || strncmp(L,'ITEM:',5)
                fseek(fid,pos,'bof');
                break;
            end
            vals = sscanf(L,'%f %f %f %f %f');
            if numel(vals)==5
                B(end+1,1:5)=vals(:).'; %#ok<SAGROW>
            end
        end
        done = true;
    end
end
fclose(fid);

assert(~isempty(B),'No numeric entries read from timestep 0 in %s', bondsDumpFile);
fprintf('Loaded %d bonds from timestep 0.\n', size(B,1));

bond_type = B(:,1);
batom1_d  = B(:,2);
batom2_d  = B(:,3);
r_d       = B(:,4);
f_d       = B(:,5);


% Build keys min(i,j)_max(i,j)
ij_min_d = min(batom1_d, batom2_d);
ij_max_d = max(batom1_d, batom2_d);

% Create string keys for matching (robust to order)
keys_dump = cell(size(ij_min_d));
for k = 1:numel(ij_min_d)
    keys_dump{k} = sprintf('%d_%d', ij_min_d(k), ij_max_d(k));
end

%% ----------- Parse bond.table -----------
fprintf('Reading %s ...\n', bondTableFile);
fid = fopen(bondTableFile,'r');
assert(fid>0, 'Could not open %s', bondTableFile);

% We'll parse line-by-line to skip comments/headers robustly
idx_t = [];  i_t = [];  j_t = [];  N_t = [];  b_t = [];
ln = fgetl(fid);
while ischar(ln)
    s = strtrim(ln);
    if isempty(s) || s(1)=='#' || strncmpi(s,'KEY',3) || (numel(s)>=1 && s(1)=='N')
        ln = fgetl(fid);
        continue;
    end
    % Expect five columns: index, i, j, N, b
    vals = sscanf(s, '%f %f %f %f %f');
    if numel(vals)==5
        idx_t(end+1,1) = vals(1); %#ok<SAGROW>
        i_t(end+1,1)   = vals(2);
        j_t(end+1,1)   = vals(3);
        N_t(end+1,1)   = vals(4);
        b_t(end+1,1)   = vals(5);
    end
    ln = fgetl(fid);
end
fclose(fid);

assert(~isempty(idx_t), 'No numeric entries were read from %s.', bondTableFile);

ij_min_t = min(i_t, j_t);
ij_max_t = max(i_t, j_t);
keys_tab = cell(size(ij_min_t));
for k = 1:numel(ij_min_t)
    keys_tab{k} = sprintf('%d_%d', ij_min_t(k), ij_max_t(k));
end

% Map from key -> [N, b]
% (containers.Map needs unique keys; if duplicates exist, last one wins)
M = containers.Map('KeyType','char','ValueType','any');
for k = 1:numel(keys_tab)
    M(keys_tab{k}) = [N_t(k), b_t(k)];
end

%% ----------- Align and compute pre-stretch -----------
nB = numel(keys_dump);
N_matched = nan(nB,1);
b_matched = nan(nB,1);

nMiss = 0;
for k = 1:nB
    key = keys_dump{k};
    if M.isKey(key)
        val = M(key);
        N_matched(k) = val(1);
        b_matched(k) = val(2);
    else
        nMiss = nMiss + 1;
    end
end

if nMiss>0
    warning('Could not find %d / %d bonds from dump in bond.table.', nMiss, nB);
end

% Keep only bonds with matching N and b
keep = ~isnan(N_matched) & ~isnan(b_matched);
r     = r_d(keep);
Nkuhn = N_matched(keep);
bkuhn = b_matched(keep);
btp   = bond_type(keep);

lambda0 = r ./ (Nkuhn .* bkuhn);

%% ----------- Quick summary -----------
fprintf('\nMatched bonds: %d / %d\n', nnz(keep), nB);
fprintf('End-to-end length r:   mean = %.3f, std = %.3f, min = %.3f, max = %.3f\n', ...
    mean(r), std(r), min(r), max(r));
fprintf('Kuhn segments N:       mean = %.2f, std = %.2f, min = %.0f, max = %.0f\n', ...
    mean(Nkuhn), std(Nkuhn), min(Nkuhn), max(Nkuhn));
fprintf('Pre-stretch lambda0:   mean = %.3f, std = %.3f, min = %.3f, max = %.3f\n\n', ...
    mean(lambda0), std(lambda0), min(lambda0), max(lambda0));

%% ----------- Plots -----------
% Common color you chose
col = [86 184 112]/255;

% 1) End-to-end length distribution (r)
figure('Name','End-to-end length distribution'); hold on;
nbins_r = max(30, min(50, round(sqrt(numel(r)))));
switch plotMode
    case 'hist'
        histogram(r, nbins_r, 'FaceColor',col, 'FaceAlpha',0.9, ...
            'EdgeColor',[0 0 0], 'Normalization','probability');
    case 'line'
        [xc_r, pdf_r] = binned_pdf_cont(r, nbins_r); % PDF using same binning spirit
        plot(xc_r, pdf_r, 'LineWidth',1.5, 'Color',col);
end
xlabel('End-to-end length r'); ylabel('Count'); % label kept per your settings
title('End-to-end Length Distribution (all bonds)');
box on; grid off;

% 2) Contour-length distribution (Kuhn segments N)
figure('Name','Contour length (Kuhn segments)'); hold on;
Nmin = floor(min(Nkuhn)); Nmax = ceil(max(Nkuhn));
edgesN = (Nmin-0.5):(Nmax+0.5); % integer-centered bins
switch plotMode
    case 'hist'
        histogram(Nkuhn, edgesN, 'FaceColor',col, 'FaceAlpha',0.85, ...
            'EdgeColor',[0 0 0], 'Normalization','probability');
    case 'line'
        [xc_N, pmf_N] = binned_pmf_discrete(Nkuhn, edgesN); % PMF at integer centers
        plot(xc_N, pmf_N, 'LineWidth',1.5, 'Color',col, 'LineStyle','-');
end
xlabel('Kuhn segments N per bond'); ylabel('Count');
title('Contour-length (N) Distribution');
box on; grid off;

% 3) Pre-stretch at equilibrium lambda0 = r/(N*b)
figure('Name','Pre-stretch distribution'); hold on;
nbins_l = max(30, min(50, round(sqrt(numel(lambda0)))));
switch plotMode
    case 'hist'
        histogram(lambda0, nbins_l, 'FaceColor',col, 'FaceAlpha',0.85, ...
            'EdgeColor',[0 0 0], 'Normalization','probability');
    case 'line'
        [xc_l, pdf_l] = binned_pdf_cont(lambda0, nbins_l);
        plot(xc_l, pdf_l, 'LineWidth',1.5, 'Color',col);
end
xlabel('\lambda_0 = r / (N \cdot b)'); ylabel('Count');
title('Pre-stretch at Equilibrium');
box on; grid off;


%% ----------- Optional: per-type overlays (uncomment to use) -----------
%{
types = unique(btp);
figure('Name','End-to-end by bond type'); hold on;
for t = reshape(types,1,[])
    rr = r(btp==t);
    histogram(rr, max(15, round(sqrt(numel(rr)))), 'DisplayStyle','stairs', 'LineWidth',1.2);
end
xlabel('r'); ylabel('Count'); title('End-to-end by bond type'); legend(arrayfun(@(x)sprintf('type %d',x), types,'uni',0));

figure('Name','Pre-stretch by bond type'); hold on;
for t = reshape(types,1,[])
    ll = lambda0(btp==t);
    histogram(ll, max(15, round(sqrt(numel(ll)))), 'DisplayStyle','stairs', 'LineWidth',1.2);
end
xlabel('\lambda_0'); ylabel('Count'); title('Pre-stretch by bond type'); legend(arrayfun(@(x)sprintf('type %d',x), types,'uni',0));
%}



