function [Nvec] = NetworkGenAssignKuhnBimodal(Bonds, options)
% NetworkGenAssignKuhnBimodal - Assign Kuhn segment numbers to bonds in a bimodal network
%
% INPUT
%   Bonds   : [BondCount x 4] -> [bondID, id1, id2, L] 
%   options : struct with fields:
%       .b
%       .bimodal.N1
%       .bimodal.N2
%
% OUTPUT
%   Nvec    : [BondCount x 1] integer vector of N per bond (same ordering as Bonds)


% ---------- Early exits / inputs ----------
Bond_count = size(Bonds,1);
Nvec = zeros(Bond_count,1);

bd = options.bimodal;
b  = options.b;

if Bond_count == 0
    return;
end

Lvec = Bonds(:,4);
type = Bonds(:,5);

mode  = options.bimodal.distribution_assignment_mode; %'single'; % 'single' | 'geom'
options.bimodal

switch lower(mode)
    case 'geom'
        % N ~ (L/b)^2 with rounding policy
        raw = (Lvec./b).^2;
        switch lower(bd.kuhn_rounding)
            case 'ceil',  Nvec = ceil(raw);
            case 'floor', Nvec = floor(raw);
            otherwise,    Nvec = round(raw);
        end
        Nvec = max(Nvec, options.bimodal.min_N);

    case 'single'
        N1 = bd.N1;
        N2 = bd.N2;
        
        Nvec(type == 1) = N1;
        Nvec(type == 2) = N2;
        
    otherwise
        error('Unknown distribution_assignment_mode: %s', mode);
end

end