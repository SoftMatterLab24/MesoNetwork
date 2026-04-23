function [Atoms, LatticeData] = AddAtoms(obj)
% -------------------------------------------------------------------------
% AddAtoms
% - Generate/scatter atoms based on obj.architecture.geometry
% - Reads all needed settings from obj.domain and obj.architecture
%
% INPUT:
%   obj : network class object
%
% OUTPUT:
%   Atoms       : atom array
%   LatticeData : [] for random geometry, struct for hex geometry
% -------------------------------------------------------------------------

    geom = lower(obj.architecture.geometry);

    switch geom

        case 'random'
            Atoms = AddAtomsRandom(obj);
            LatticeData = [];

        case 'hex_lattice'
            [Atoms, LatticeData] = AddAtomsHex(obj);

        case 'bottle_brush'
            Atoms = AddAtomsBottleBrush(obj);
            LatticeData = [];

        otherwise
            error('AddAtoms: unknown geometry "%s".', obj.architecture.geometry);

    end

end