function syncKuhnAssignmentFromTopology(obj)
% -------------------------------------------------------------------------
% syncKuhnAssignmentFromTopology
% - If auto mode is set to auto, inherit the strand-topology
%   assignment settings from obj.architecture.strand_typology.
% - Copies only the active mode block and sets obj.perbond.kuhn.mode to the
%   same distribution type.
%
% INPUT:
%   obj : network object
%
% OUTPUT:
%   none
%   (obj.perbond.kuhn is modified in place)
% -------------------------------------------------------------------------

    if ~obj.architecture.strand_typology.auto
        return;
    end

    topo = obj.architecture.strand_typology;
    mode = lower(topo.mode);

    obj.perbond.kuhn.mode = mode;

    switch mode
        case 'mono'
            obj.perbond.kuhn.mono = topo.mono;

        case 'uniform'
            obj.perbond.kuhn.uniform = topo.uniform;

        case {'poly', 'polydisperse'}
            obj.perbond.kuhn.poly = topo.poly;

        case 'bimodal'
            obj.perbond.kuhn.bimodal = topo.bimodal;

        otherwise
            error('syncKuhnAssignmentFromTopology: unknown strand typology mode "%s".', topo.mode);
    end
end