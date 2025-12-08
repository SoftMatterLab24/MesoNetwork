function Nvec = NetworkGenAssignKuhnBimodal(Bonds, Atoms, options)
% NetworkGenAssignKuhnBimodal - Assign Kuhn segment numbers to bonds in a bimodal network
%
% INPUT
%   Bonds   : [BondCount x 5] -> [bondID, id1, id2, L, type]
%   options :
%       .b
%       .bimodal.N1
%       .bimodal.N2
%       .bimodal.min_N                 (default 1)
%       .bimodal.N_cap                 (optional; upper cap)
%       .bimodal.distribution_assignment_mode  ('single'|'geom'|'gaussian')
%       .bimodal.kuhn_rounding         ('round'|'ceil'|'floor')   % used in 'geom'
%        % For 'typed_gaussian' you may specify widths by std OR FWHM:
%       .bimodal.gauss_std1,  .bimodal.gauss_std2
%       .bimodal.gauss_fwhm1, .bimodal.gauss_fwhm2
%
% OUTPUT
%   Nvec  : [BondCount x 1] integer N per bond
natoms      = size(Atoms,1);
Bond_count = size(Bonds,1);
Nvec = zeros(Bond_count,1);
if Bond_count == 0, return; end

bd = options.bimodal;
b  = options.b;

Lvec = Bonds(:,4);
type = Bonds(:,5);

% defaults
if ~isfield(bd,'min_N'),   bd.min_N = 1; end
has_cap = isfield(bd,'N_cap') && ~isempty(bd.N_cap);
is_manual = strcmpi(bd.bin_window_method,'manual');

mode = lower(bd.distribution_assignment_mode);

switch mode
    case 'geom'
        % N ~ (L/b)^2 with rounding policy
        raw = (Lvec./max(b,eps)).^2;
        if isfield(bd,'kuhn_rounding')
            switch lower(bd.kuhn_rounding)
                case 'ceil',  Nvec = ceil(raw);
                case 'floor', Nvec = floor(raw);
                otherwise,    Nvec = round(raw);
            end
        else
            Nvec = round(raw);
        end
        Nvec = max(Nvec, bd.min_N);
        if has_cap, Nvec = min(Nvec, bd.N_cap); end

    case 'gaussian'
        % assign N by Gaussian distribution with long chains with larger N
        N1 = bd.N1;  N2 = bd.N2; % Per-type Gaussian around N1/N2 with small widths

        % Resolve widths (std priority, else FWHM/2.355, else sensible defaults)
        if isfield(bd,'std1') && ~isempty(bd.std1) && is_manual
            s1 = bd.std1;
        elseif isfield(bd,'gauss_fwhm1') && ~isempty(bd.gauss_fwhm1)
            s1 = bd.gauss_fwhm1 / 2.355;
        else
            s1 = max(1, 0.15*N1);   % default ~15% of mean
        end

        if isfield(bd,'std2') && ~isempty(bd.std2) && is_manual
            s2 = bd.std2;
        elseif isfield(bd,'gauss_fwhm2') && ~isempty(bd.gauss_fwhm2)
            s2 = bd.gauss_fwhm2 / 2.355;
        else
            s2 = max(1, 0.20*N2);   % default ~20% of mean
        end

        % Sample typed normals, integerize, and clamp
        idx1 = (type == 1); n1 = sum(idx1);
        idx2 = (type == 2); n2 = sum(idx2);

        if n1 > 0
            n_draw1 = N1 + s1 * randn(n1,1);
            n_draw1 = round(n_draw1);
            n_draw1 = max(n_draw1, bd.min_N);
            if has_cap, n_draw1 = min(n_draw1, bd.N_cap); end

            N_list = sort(n_draw1,'ascend'); 
            
            [~, idxL] = sort(Lvec(idx1), 'ascend');
            
            % assign largest N to longest bonds
            n_draw1_sorted(idxL) = N_list;
            Nvec(idx1) = n_draw1_sorted;
        
        end
        if n2 > 0
            n_draw2 = N2 + s2 * randn(n2,1);
            n_draw2 = round(n_draw2);
            n_draw2 = max(n_draw2, bd.min_N);
            if has_cap, n_draw2 = min(n_draw2, bd.N_cap); end
            Nvec(idx2) = n_draw2;

            N_list = sort(n_draw2,'ascend');
            [~, idxL] = sort(Lvec(idx2), 'ascend'); 

            n_draw2_sorted(idxL) = N_list;
            Nvec(idx2) = n_draw2_sorted;
        end

        % Edge case: if any zeros slipped in (e.g., NaN), repair
        bad = ~isfinite(Nvec) | (Nvec < bd.min_N);
        if any(bad)
            Nvec(bad) = bd.min_N;
        end

        % Check end-to-end vs contour length
        notEnoughSegments = Lvec./(Nvec*b) > 1;
        % Overide: Use geometric scaling for chains that dont have enough
        if any(notEnoughSegments)
            Nvec(notEnoughSegments) = ceil( 1.1 * (Lvec(notEnoughSegments)/b) );
        end
        
    case 'single'
        % Assign exactly N1/N2 by stored type
        N1 = bd.N1;  N2 = bd.N2;
        Nvec(type == 1) = N1;
        Nvec(type == 2) = N2;
        Nvec = max(Nvec, bd.min_N);
        if has_cap, Nvec = min(Nvec, bd.N_cap); end
    otherwise
        error('Unknown distribution_assignment_mode: %s', bd.distribution_assignment_mode);
    end

    fprintf("Kuhn-to-crosslinker ratio %0.4f \n",sum(Nvec)/natoms);
    fprintf("Average chain length %0.4f \n",sum(Nvec)/length(Nvec))
end
