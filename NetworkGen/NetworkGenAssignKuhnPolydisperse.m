function Nvec = NetworkGenAssignKuhnPolydisperse(Bonds, options)
%NETWORKGENASSIGNKUHNPOLYDISPERSE
% Assign per-bond Kuhn segment counts N for a polydisperse network.
% 
% INPUT
%   Bonds   : [BondCount x 4] -> [bondID, id1, id2, L]
%   options : struct with fields:
%       .b
%       .polydisperse.distribution_assignment_mode  ('geom'|'range'|'pmf')
%       .polydisperse.min_N
%       .polydisperse.align_to_length               ('ascend'|'none')
%       .polydisperse.kuhn_rounding                 ('round'|'ceil'|'floor')
%       .polydisperse.plot_hist                     (true|false)
%       .polydisperse.N_range_method                ('rank'|'linear')
%       .polydisperse.N_target_min                  (int)
%       .polydisperse.N_target_max                  (int)
%       .polydisperse.pmf_nu0
%       .polydisperse.pmf_meanN
%       .polydisperse.pmf_cut_mode                  ('cap')
%       .polydisperse.pmf_nu_max
%       .polydisperse.integerize_rule               ('largest_remainder')
%
% OUTPUT
%   Nvec    : [BondCount x 1] integer vector of N per bond (same ordering as Bonds)

% ---------- Early exits / inputs ----------
Bond_count = size(Bonds,1);
Nvec = zeros(Bond_count,1);
if Bond_count == 0
    return;
end

pd = options.polydisperse;
b  = options.b;

mode   = pd.distribution_assignment_mode;
min_N  = pd.min_N;

Lvec = Bonds(:,4);

% ---------- Mode branches ----------
switch lower(mode)
    case 'geom'
        % N ~ (L/b)^2 with rounding policy
        raw = (Lvec./b).^2;
        switch lower(pd.kuhn_rounding)
            case 'ceil',  Nvec = ceil(raw);
            case 'floor', Nvec = floor(raw);
            otherwise,    Nvec = round(raw);
        end
%         Nvec = max(Nvec, min_N);
        Nvec = 

    case 'range'
        % Map lengths monotonically to [N_target_min, N_target_max]
        Nlo = min(pd.N_target_min, pd.N_target_max);
        Nhi = max(pd.N_target_min, pd.N_target_max);

        meth = lower(pd.N_range_method);
        if strcmp(meth,'rank')
            [Ls, idx] = sort(Lvec, 'ascend'); %#ok<ASGLU>
            if Bond_count == 1
                Ntargets = (Nlo + Nhi)/2;
            else
                Ntargets = linspace(Nlo, Nhi, Bond_count).';
            end
            Ntmp = zeros(Bond_count,1);
            Ntmp(idx) = Ntargets;
            Nvec = round(Ntmp);
        else
            Lmin = min(Lvec); Lmax = max(Lvec);
            if Lmax == Lmin
                Nvec = round(((Nlo + Nhi)/2) * ones(Bond_count,1));
            else
                t = (Lvec - Lmin) ./ (Lmax - Lmin);
                Nvec = round(Nlo + t .* (Nhi - Nlo));
            end
        end
        Nvec = max(Nvec, min_N);
        Nvec = min(Nvec, max(pd.N_target_min, pd.N_target_max));

    case 'pmf'
        % Truncated geometric on nu ∈ [nu0, nuMax], PMF �? p (1-p)^(nu-nu0)
        nu0   = pd.pmf_nu0;
        nuMax = max(pd.pmf_nu0, pd.pmf_nu_max); % ensure cap ≥ base
        K     = nuMax - nu0;                    % support in k = 0..K
        targetMeanN = pd.pmf_meanN;
        targetMeanK = max(0, targetMeanN - nu0);

        % Solve for p in (0,1) so truncated mean(k) matches targetMeanK
        % mean_k(p) = [ (1-p)*(1 - (K+1)(1-p)^K + K(1-p)^(K+1)) / p ] / [ 1 - (1-p)^(K+1) ]
        f = @(p) mean_k_of_p(p, K) - targetMeanK;

        p_lo = 1e-8; p_hi = 1-1e-8;
        Ps = linspace(1e-6, 1-1e-6, 200);
        Fs = zeros(size(Ps));
        for ii = 1:numel(Ps), Fs(ii) = f(Ps(ii)); end

        bracket_found = false; a = NaN; b = NaN;
        for ii = 1:(numel(Ps)-1)
            if Fs(ii) == 0
                a = Ps(ii); b = Ps(ii); bracket_found = true; break;
            elseif Fs(ii)*Fs(ii+1) < 0
                a = Ps(ii); b = Ps(ii+1); bracket_found = true; break;
            end
        end

        if bracket_found
            if a == b
                p_opt = a;
            else
                p_opt = fzero(f, [a, b]);
            end
        else
            % fallback: uncapped geometric mean k = (1-p)/p  =>  p ≈ 1/(targetMeanK+1)
            a_guess = max(eps, targetMeanK);
            p_opt   = 1/(a_guess + 1);
            p_opt   = min(max(p_opt, p_lo), p_hi);
        end

        p = min(max(p_opt, p_lo), p_hi);
        r = 1 - p;

        % normalized PMF on k=0..K: Pk = [p r^k] / [1 - r^(K+1)]
        denom = 1 - r^(K+1);
        Pk = (p * (r .^ (0:K)).') / denom; % column vector

        % Convert expected counts -> integer counts (largest remainder)
        exp_counts = Pk * Bond_count;
        base_counts = floor(exp_counts);
        remainder   = exp_counts - base_counts;
        assigned = sum(base_counts);
        deficit  = Bond_count - assigned;

        if deficit > 0
            [~, order] = sort(remainder, 'descend');
            for t = 1:deficit
                base_counts(order(t)) = base_counts(order(t)) + 1;
            end
        elseif deficit < 0
            [~, order] = sort(remainder, 'ascend');
            for t = 1:(-deficit)
                idx = order(t);
                if base_counts(idx) > 0
                    base_counts(idx) = base_counts(idx) - 1;
                end
            end
        end

        % Build nondecreasing ν-list according to counts
        N_list = zeros(Bond_count,1);
        ptr = 1;
        for k = 0:K
            cnt = base_counts(k+1);
            if cnt <= 0, continue; end
            val = nu0 + k;
            N_list(ptr:ptr+cnt-1) = val;
            ptr = ptr + cnt;
        end
        if ptr <= Bond_count
            N_list(ptr:Bond_count) = nu0 + K; % safety
        end

        % Align to geometry if requested
        align = 'ascend';
        if isfield(pd,'align_to_length'), align = lower(pd.align_to_length); end
        if strcmp(align,'ascend')
            [~, idxL] = sort(Lvec, 'ascend');
            Nvec(idxL) = N_list;
        else % 'none' -> random permutation to avoid bias
            rp = randperm(Bond_count).';
            Nvec(rp) = N_list;
        end

        % guard
        Nvec = max(Nvec, min_N);

    otherwise
        error('Unknown distribution_assignment_mode: %s', mode);
end

end

% ===== helper (nested at EOF for R2016a) =====
function mk = mean_k_of_p(p, K)
% Truncated geometric on k=0..K with Pk �? p (1-p)^k.
% Returns mean(k).
r = 1 - p;
num = (1-p) .* (1 - (K+1)*r.^K + K*r.^(K+1)) ./ p;
den = 1 - r.^(K+1);
mk = num ./ den;
end
