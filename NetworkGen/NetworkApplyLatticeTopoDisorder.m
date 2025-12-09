function [Atoms, Bonds] = NetworkApplyLatticeTopoDisorder(Atoms, Bonds, options)
% -------------------------------------------------------------------------
% NetworkApplyLatticeTopoDisorder
%
% Introduce topological defects by randomly deleting bonds with a
% probability that scales smoothly with:
%
%   options.lattice.disorder_level        in [0,1]
%   options.lattice.max_topo_del_per_node (e.g. 2)
%   options.lattice.min_degree_keep       (e.g. 2)
%
% Logic:
%   - Compute deletion probability per bond:
%         p_del ~ disorder_level * max_topo_del_per_node / 6
%     (6 is the typical degree in a triangular lattice.)
%   - Loop over bonds in random order, and for each bond:
%         if rand < p_del AND both endpoints have degree > min_degree_keep
%         -> delete the bond.
%
% At disorder_level = 0:
%   -> p_del = 0, no bonds removed.
%
% As disorder_level increases towards 1:
%   -> expected deletions per node increase smoothly, not in jumps.
%
% Bonds is (Nbonds x 5): [id  i  j  L0  type]
% Atoms is passed through unchanged.
% -------------------------------------------------------------------------

if ~isfield(options,'lattice')
    return;
end

if ~isfield(options.lattice,'disorder_level')
    return;
end

disorder_level = options.lattice.disorder_level;
if disorder_level <= 0
    % No topological changes
    return;
end

% Switch must be enabled
if ~isfield(options.lattice,'enable_topo_disorder')
    return;
end
if ~options.lattice.enable_topo_disorder
    return;
end

% Parameters
max_del_per_node = 2;
if isfield(options.lattice,'max_topo_del_per_node')
    max_del_per_node = options.lattice.max_topo_del_per_node;
end

min_degree_keep = 2;
if isfield(options.lattice,'min_degree_keep')
    min_degree_keep = options.lattice.min_degree_keep;
end

if max_del_per_node <= 0
    return;
end

Natoms = size(Atoms,1);
Nbonds = size(Bonds,1);
if Nbonds == 0
    return;
end

% Endpoints
i_list = Bonds(:,2);
j_list = Bonds(:,3);

% Current degree of each node
deg = zeros(Natoms,1);
k = 1;
while k <= Nbonds
    ii = i_list(k);
    jj = j_list(k);
    if ii >= 1 && ii <= Natoms
        deg(ii) = deg(ii) + 1;
    end
    if jj >= 1 && jj <= Natoms
        deg(jj) = deg(jj) + 1;
    end
    k = k + 1;
end

% -------------------------------------------------------------------------
% Smooth mapping: expected deletions per node ~ disorder_level * max_del_per_node
% For a triangular lattice, degree ~ 6, bonds ~ 3 per node.
% p_del is per-bond deletion probability:
%   total deletions ~ p_del * Nbonds
%   expected deletions per node ~ 2 * total_deletions / Natoms
%   => p_del ~ disorder_level * max_del_per_node / 6
% -------------------------------------------------------------------------
lambda_node = disorder_level * max_del_per_node;
p_del = lambda_node / 6.0;

% Clamp to [0, 0.9] to avoid pathological behavior
if p_del < 0
    p_del = 0;
elseif p_del > 0.9
    p_del = 0.9;
end

% Random visiting order over bonds
perm = randperm(Nbonds);
delete_flag = false(Nbonds,1);

idx = 1;
while idx <= Nbonds
    bID = perm(idx);
    
    if delete_flag(bID)
        idx = idx + 1;
        continue;
    end
    
    ii = i_list(bID);
    jj = j_list(bID);
    
    % Check degree constraint first
    if (deg(ii) <= min_degree_keep) || (deg(jj) <= min_degree_keep)
        idx = idx + 1;
        continue;
    end
    
    % Probabilistic deletion
    if rand < p_del
        delete_flag(bID) = true;
        deg(ii) = deg(ii) - 1;
        deg(jj) = deg(jj) - 1;
    end
    
    idx = idx + 1;
end

% Apply deletions
keep_mask = ~delete_flag;
Bonds = Bonds(keep_mask,:);

% Renumber bond IDs sequentially
Nbonds_new = size(Bonds,1);
if Nbonds_new > 0
    Bonds(:,1) = (1:Nbonds_new).';
end

end
