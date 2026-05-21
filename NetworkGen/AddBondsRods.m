function RodBonds = AddBondsRods(obj, Atoms)
% -------------------------------------------------------------------------
% AddBondsRods
% - Build sequential intra-rod bonds for bottle-brush architecture.
% - Rod membership is inferred from the molID column (Atoms(:,2)): all atoms
%   sharing a molID are one rod, bonded in the order they appear in Atoms.
% - Bond length L0 = obj.architecture.bottlebrush.sigma_r for every rod bond.
% - Bond types are NOT assigned here. The 5th column is set to 0 so that
%   AssignMultiType can label rod vs crosslink bonds by endpoint molID match.
%
% INPUT:
%   obj   : network object
%   Atoms : atom array (ID | molID | X | Y | Z | deg | nbrs...)
%
% OUTPUT:
%   RodBonds : [bondID | atomID_i | atomID_j | L0 | type=0]
%              bondID is assigned within RodBonds (1..Nrb); the AddBonds
%              dispatcher renumbers after concatenating with crosslink bonds.
% -------------------------------------------------------------------------

    sigma_r = obj.architecture.bottlebrush.sigma_r;

    ids       = Atoms(:, 1);
    mol_ids   = Atoms(:, 2);
    natom     = size(Atoms, 1);

    if natom < 2
        RodBonds = zeros(0, 5);
        return;
    end

    % Allocate with upper bound. Worst case: every atom sequential with its
    % neighbor -> natom-1 bonds (overestimate; actual is natom - N_rods).
    RodBonds = zeros(natom, 5);
    nbond = 0;

    % Build rod bonds by grouping atoms by molID and connecting consecutive
    % rod atoms in ascending ID order. Rod atoms are created with incremental
    % IDs, so this prevents same-molID bonds from spanning across nonadjacent
    % rod atoms.
    unique_mol_ids = unique(mol_ids(mol_ids > 0));
    for mi = 1:numel(unique_mol_ids)
        mol = unique_mol_ids(mi);
        idx = find(mol_ids == mol);
        if numel(idx) < 2
            continue;
        end

        % Order rod atoms by atom ID rather than by nearest spatial neighbor.
        [~, order] = sort(ids(idx));
        chain = idx(order);

        for j = 1:(numel(chain) - 1)
            nbond = nbond + 1;
            RodBonds(nbond, :) = [nbond, ids(chain(j)), ids(chain(j+1)), sigma_r, 0];
        end
    end

    RodBonds = RodBonds(1:nbond, :);

    obj.log.print('   AddBondsRods: built %d intra-rod bonds across %d unique molIDs\n', ...
        nbond, numel(unique(mol_ids(mol_ids > 0))));

end
