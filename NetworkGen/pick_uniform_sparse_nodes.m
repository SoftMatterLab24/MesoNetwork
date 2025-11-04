function isSparse = pick_uniform_sparse_nodes(x, y, f_sparse, r_spacing)
    % pick_uniform_sparse_nodes: selects uniformly spaced random subset of nodes
    %  x,y         - node coordinates
    %  f_sparse    - fraction of nodes desired (e.g. 0.2)
    %  r_spacing   - minimum spacing between selected nodes (in same units as x,y)

    natom = numel(x);
    Nsparse_target = round(f_sparse * natom);
    isSparse = false(natom,1);

    idx_all = randperm(natom);
    picked = [];

    for k = 1:natom
        i = idx_all(k);
        if isempty(picked)
            picked(end+1) = i;
            continue;
        end
        % Check distance to all previously selected
        dx = x(picked) - x(i);
        dy = y(picked) - y(i);
        if all( sqrt(dx.^2 + dy.^2) >= r_spacing )
            picked(end+1) = i;
        end
        if numel(picked) >= Nsparse_target
            break;
        end
    end
    isSparse(picked) = true;
end
