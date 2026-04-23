function Nvec = AssignPerBond(obj, Bonds, Atoms)
% -------------------------------------------------------------------------
% AssignPerBond
% - Top-level dispatcher for per-bond Kuhn segment assignment.
%
% This version is type-agnostic: every bond is processed by the same
% distribution implied by obj.perbond.kuhn.mode, regardless of Bonds(:,5).
% Multi-type per-bond assignment (e.g. different distributions for
% different bond types) is deferred to a future version.
%
% INPUT:
%   obj   : network object
%   Bonds : bond array [bondID, id1, id2, L0, type]
%   Atoms : atom array
%
% OUTPUT:
%   Nvec  : [Nbonds x 1] Kuhn segment count per bond
% -------------------------------------------------------------------------

    nbonds = size(Bonds, 1);

    if nbonds == 0
        Nvec = zeros(0, 1);
        return;
    end

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
