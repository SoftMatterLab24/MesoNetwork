function [Atoms, Bonds] = AddBonds(obj, Atoms, LatticeData)
% -------------------------------------------------------------------------
% AddBonds
% - Top-level dispatcher for bond creation
% - Decides internal bond-generation routine from:
%       obj.architecture.strand_typology.mode
%       obj.architecture.geometry
%
% INPUT:
%   obj         : network object
%   Atoms       : atom array
%   LatticeData : [] for random geometry, struct for lattice geometry
%
% OUTPUT:
%   Atoms : updated atom array with degree / neighbor list rebuilt
%   Bonds : bond array
% -------------------------------------------------------------------------

    mode = lower(obj.architecture.strand_typology.mode);

    switch mode

        case 'mono'
            [Atoms, Bonds] = AddBondsMono(obj, Atoms, LatticeData);

        case {'poly', 'polydisperse'}
            [Atoms, Bonds] = AddBondsPoly(obj, Atoms, LatticeData);

        case 'bimodal'
            [Atoms, Bonds] = AddBondsBimodal(obj, Atoms, LatticeData);

        otherwise
            error('AddBonds: unknown strand_typology.mode "%s".', ...
                  obj.architecture.strand_typology.mode);
    end

end