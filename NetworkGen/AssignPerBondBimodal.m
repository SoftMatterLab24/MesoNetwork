function Nvec = AssignPerBondBimodal(obj, Bonds, Atoms)
% -------------------------------------------------------------------------
% AssignPerBondBimodal
% - Assign Kuhn segment numbers to bonds in a bimodal manner.
%
% Per the type-agnostic rule for AssignPerBond, the bond-type column is
% IGNORED. Bonds are split into two groups by median length:
%       shorter half -> group 1 (target mean N1)
%       longer half  -> group 2 (target mean N2)
%
% Supported methods (obj.perbond.kuhn.bimodal.method):
%       'single'    -> constant N1 for group 1, N2 for group 2
%       'geom'      -> N ~ (L/b)^2 per bond, rounded
%       'gaussian'  -> typed Gaussian around N1 / N2 with length-sorted
%                      placement within each group
%
% OUTPUT:
%   Nvec : [Nbonds x 1]
% -------------------------------------------------------------------------

    natoms = size(Atoms,1);
    nbonds = size(Bonds,1);

    Nvec = zeros(nbonds,1);
    if nbonds == 0
        return;
    end

    bd = obj.perbond.kuhn.bimodal;
    b  = obj.domain.b;

    Lvec = Bonds(:,4);

    % Length-based split (bond-type column ignored by design)
    type = ones(nbonds,1);
    if nbonds > 1
        Lmid = median(Lvec);
        type(Lvec > Lmid) = 2;
    end

    if isempty(bd.mean_1)
        bd.mean_1 = 35;
    end
    if isempty(bd.mean_2)
        bd.mean_2 = 60;
    end

    N1 = round(bd.mean_1);
    N2 = round(bd.mean_2);

    if isfield(bd,'min_value') && ~isempty(bd.min_value)
        min_N = round(bd.min_value);
    else
        min_N = 1;
    end
    min_N = max(min_N,1);

    has_cap = false;
    N_cap = [];
    if isfield(bd,'target_max') && ~isempty(bd.target_max)
        has_cap = true;
        N_cap = round(bd.target_max);
    end

    mode = lower(bd.method);

    switch mode

        case 'geom'
            raw = (Lvec ./ max(b,eps)).^2;

            if isfield(bd,'rounding') && ~isempty(bd.rounding)
                switch lower(bd.rounding)
                    case 'ceil'
                        Nvec = ceil(raw);
                    case 'floor'
                        Nvec = floor(raw);
                    otherwise
                        Nvec = round(raw);
                end
            else
                Nvec = round(raw);
            end

            Nvec = max(Nvec, min_N);
            if has_cap
                Nvec = min(Nvec, N_cap);
            end

        case 'gaussian'
            is_manual = false;
            if isfield(bd,'bin_window_method') && ~isempty(bd.bin_window_method)
                is_manual = strcmpi(bd.bin_window_method,'manual');
            end

            if isfield(bd,'std_1') && ~isempty(bd.std_1) && is_manual
                s1 = bd.std_1;
            else
                s1 = max(1, 0.15*N1);
            end

            if isfield(bd,'std_2') && ~isempty(bd.std_2) && is_manual
                s2 = bd.std_2;
            else
                s2 = max(1, 0.20*N2);
            end

            idx1 = (type == 1);
            idx2 = (type == 2);

            n1 = sum(idx1);
            n2 = sum(idx2);

            if n1 > 0
                draw1 = N1 + s1 * randn(n1,1);
                draw1 = round(draw1);
                draw1 = max(draw1, min_N);
                if has_cap
                    draw1 = min(draw1, N_cap);
                end

                N_list = sort(draw1, 'ascend');
                [~, idxL] = sort(Lvec(idx1), 'ascend');

                tmp = zeros(n1,1);
                tmp(idxL) = N_list;
                Nvec(idx1) = tmp;
            end

            if n2 > 0
                draw2 = N2 + s2 * randn(n2,1);
                draw2 = round(draw2);
                draw2 = max(draw2, min_N);
                if has_cap
                    draw2 = min(draw2, N_cap);
                end

                N_list = sort(draw2, 'ascend');
                [~, idxL] = sort(Lvec(idx2), 'ascend');

                tmp = zeros(n2,1);
                tmp(idxL) = N_list;
                Nvec(idx2) = tmp;
            end

            bad = ~isfinite(Nvec) | (Nvec < min_N);
            if any(bad)
                Nvec(bad) = min_N;
            end

            % Ensure contour length can exceed end-to-end distance
            notEnoughSegments = (Lvec ./ (Nvec * b)) > 1;
            if any(notEnoughSegments)
                Nvec(notEnoughSegments) = ceil(1.1 * (Lvec(notEnoughSegments) / b));
            end

        case 'single'
            Nvec(type == 1) = N1;
            Nvec(type == 2) = N2;
            Nvec = max(Nvec, min_N);
            if has_cap
                Nvec = min(Nvec, N_cap);
            end

        otherwise
            error('AssignPerBondBimodal: unknown bimodal.method "%s".', bd.method);
    end

    obj.log.print('   Kuhn assignment mode: bimodal (%s, length-split)\n', mode);
    obj.log.print('   Kuhn-to-crosslinker ratio %0.4f\n', sum(Nvec)/max(natoms,1));
    obj.log.print('   Average chain length %0.4f\n', sum(Nvec)/max(length(Nvec),1));

end
