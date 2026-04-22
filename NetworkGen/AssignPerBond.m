function Nvec = AssignPerBond(obj, Bonds, Atoms)
% -------------------------------------------------------------------------
% AssignPerBond
% - Top-level dispatcher for per-bond Kuhn segment assignment
%
% LOGIC:
%   1) If obj.perbond.kuhn.auto == true:
%         use obj.architecture.strand_typology.mode
%   2) Else:
%         use obj.perbond.kuhn.mode
%
% This allows mixed cases, e.g.
%   - mono topology + poly Kuhn assignment
%   - mono lattice + bimodal Kuhn assignment
%
% INPUT:
%   obj   : network object
%   Bonds : bond array
%   Atoms : atom array
%
% OUTPUT:
%   Nvec  : [Nbonds x 1] Kuhn segment count per bond
% -------------------------------------------------------------------------

    mode = lower(obj.perbond.kuhn.mode);

    switch mode

        case 'mono'
            Nvec = AssignPerBondMono(obj, Bonds, Atoms);

        case {'poly', 'polydisperse'}
            Nvec = AssignPerBondPoly(obj, Bonds, Atoms);

        case 'bimodal'
            Nvec = AssignPerBondBimodal(obj, Bonds, Atoms);

        otherwise
            error('AssignPerBond: unknown Kuhn assignment mode "%s".', mode);
    end

end