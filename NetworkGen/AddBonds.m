function [Atoms, Bonds] = AddBonds(obj, Atoms, LatticeData)
% -------------------------------------------------------------------------
% AddBonds
% - Top-level bond-creation dispatcher.
%
% For bottle_brush geometry: first build intra-rod bonds via AddBondsRods,
% then hand them to the chosen crosslink topology (mono / poly / bimodal)
% as PreBonds so adjacency is seeded and duplicates are avoided. The
% Max_bond cap is temporarily inflated by Nr for the crosslink phase,
% matching the old bottle-brush behavior.
%
% For other geometries: rod bonds are empty; crosslink topology runs alone.
%
% No routine assigns bond types here. Bonds(:,5) is left 0 for every bond
% created; AssignMultiType writes the real type in a later pipeline stage.
%
% INPUT:
%   obj         : network object
%   Atoms       : atom array (new layout with molID at col 2)
%   LatticeData : [] for random / bottle_brush, struct for hex_lattice
%
% OUTPUT:
%   Atoms : atom array with degree/neighbor lists rebuilt on full bond set
%   Bonds : [bondID | atomID_i | atomID_j | L0 | type=0]
% -------------------------------------------------------------------------

    geom = lower(obj.architecture.geometry);
    Max_peratom_bond = obj.peratom.Max_peratom_bond;

    % -----------------------------------------------------------------
    % Phase 1: rod bonds (bottle_brush only)
    % -----------------------------------------------------------------
    if strcmpi(geom, 'bottle_brush')
        RodBonds = AddBondsRods(obj, Atoms);

        % Inflate Max_bond for the crosslink phase so the loop doesn't
        % terminate early on the rod-inflated atom count. Restore after.
        Nr = obj.architecture.bottlebrush.Nr;
        original_Max_bond = obj.domain.Max_bond;
        obj.domain.Max_bond = Nr * original_Max_bond;
    else
        RodBonds = zeros(0, 5);
        original_Max_bond = [];
    end

    % -----------------------------------------------------------------
    % Phase 2: crosslink topology (runs for every geometry)
    % -----------------------------------------------------------------
    mode = lower(obj.architecture.strand_typology.mode);

    switch mode

        case 'mono'
            [Atoms, Bonds] = AddBondsMono(obj, Atoms, LatticeData, RodBonds);

        case {'poly', 'polydisperse'}
            [Atoms, Bonds] = AddBondsPoly(obj, Atoms, LatticeData, RodBonds);

        case 'bimodal'
            [Atoms, Bonds] = AddBondsBimodal(obj, Atoms, LatticeData, RodBonds);

        otherwise
            error('AddBonds: unknown strand_typology.mode "%s".', ...
                  obj.architecture.strand_typology.mode);
    end

    % Restore Max_bond if we mutated it
    if ~isempty(original_Max_bond)
        obj.domain.Max_bond = original_Max_bond;
    end

    % -----------------------------------------------------------------
    % Phase 3: combine rod + crosslink bonds, renumber, rebuild neighbors
    % -----------------------------------------------------------------
    if ~isempty(RodBonds)
        Bonds = [RodBonds; Bonds];
        if ~isempty(Bonds)
            Bonds(:,1) = (1:size(Bonds,1)).';
        end

        % Crosslink finalize skipped internal pruning (skip_prune=true for
        % bottle_brush via PreBonds), so we rebuild neighbor lists here on
        % the combined bond set.
        Atoms = rebuild_atom_neighbors(Atoms, Bonds, Max_peratom_bond);
    end

end


% =========================================================================
% Helper: rebuild degree / neighbor lists on the new column layout
% =========================================================================
function Atoms = rebuild_atom_neighbors(Atoms, Bonds, Max_peratom_bond)

    needed_cols = 6 + Max_peratom_bond;
    if size(Atoms,2) < needed_cols
        Atoms(:, size(Atoms,2)+1 : needed_cols) = 0;
    end

    Atoms(:,6) = 0;
    Atoms(:,7:6+Max_peratom_bond) = 0;

    for k = 1:size(Bonds,1)
        ii = Bonds(k,2);
        jj = Bonds(k,3);

        nb1 = Atoms(ii,6) + 1;
        if nb1 <= Max_peratom_bond
            Atoms(ii,6) = nb1;
            Atoms(ii, 6 + nb1) = jj;
        end

        nb2 = Atoms(jj,6) + 1;
        if nb2 <= Max_peratom_bond
            Atoms(jj,6) = nb2;
            Atoms(jj, 6 + nb2) = ii;
        end
    end

end
