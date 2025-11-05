function [Domain,Atoms,Bonds] = NetworkScaleDomain(Domain,Atoms,Bonds,scaleFactor)
% -------------------------------------------------------------------------
% NetworkScaleDomain
% Uniformly rescales the domain bounds, atom coordinates, and bond lengths.
%
% INPUT:
%   Domain: structure with fields xlo,xhi,ylo,yhi,zlo,zhi
%   Atoms:  [ID, x, y, z, degree, ...]
%   Bonds:  [ID, i, j, L, type] (optional)
%   scaleFactor: scalar multiplier (e.g., 0.5 halves domain size)
%
% OUTPUT:
%   Domain, Atoms, Bonds: rescaled versions
% -------------------------------------------------------------------------

if nargin < 4
    error('Usage: [Domain,Atoms,Bonds] = NetworkScaleDomain(Domain,Atoms,Bonds,scaleFactor)');
end

% --- Scale domain bounds ---
Domain.xlo = Domain.xlo * scaleFactor;
Domain.xhi = Domain.xhi * scaleFactor;
Domain.ylo = Domain.ylo * scaleFactor;
Domain.yhi = Domain.yhi * scaleFactor;
% Domain.zlo = Domain.zlo * scaleFactor;
% Domain.zhi = Domain.zhi * scaleFactor;

% --- Scale atom coordinates ---
Atoms(:,2:4) = Atoms(:,2:4) * scaleFactor;

% --- Scale bond lengths (4th column) if present ---
if ~isempty(Bonds) && size(Bonds,2) >= 4
    Bonds(:,4) = Bonds(:,4) * scaleFactor;
end

end
