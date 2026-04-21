function network_prefix = resolveNetworkTypePrefix(obj)
% -------------------------------------------------------------------------
% resolveNetworkTypePrefix
% - Returns:
%       'BD' for bimodal
%       'PD' for polydisperse
%       'MD' for monodisperse
%
% Logic is based on topology mode primarily.
% -------------------------------------------------------------------------

    topo_mode = lower(obj.architecture.strand_typology.mode);

    switch topo_mode
        case 'bimodal'
            network_prefix = 'BD';

        case {'poly', 'polydisperse'}
            network_prefix = 'PD';

        otherwise
            network_prefix = 'MD';
    end

end