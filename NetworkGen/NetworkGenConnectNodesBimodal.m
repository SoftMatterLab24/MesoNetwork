function [Atoms,Bonds] = NetworkGenConnectNodesBimodal(Domain,Atoms,options,advancedOptions)
% NetworkGenConnectNodesBimodal - Connect Atoms nodes with Bonds for bimodal distribution
% INPUT:
%   Atoms layout: [ ID | X | Y | Z | num_bond | nbr1 | nbr2 | nbr3 | nbr4 | ... ]
%   IDs are arbitrary (NOT equal to row index).
% OUTPUT:
%   BondsOut: [bondID | id1 | id2 | L | bond_type]  (ids, not rows)
%   AtomsOut: Atoms with num_bond and neighbor slots rebuilt to match BondsOut

% --------- Unpack ---------
natom             = size(Atoms,1);
Max_bond          = Domain.Max_bond;
Max_peratom_bond  = Domain.Max_peratom_bond;
global_limit      = Domain.bond_global_try_limit;
stall_limit       = Domain.max_attempts_without_progress;
min_keep          = Domain.min_degree_keep;

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;
zlo = Domain.zlo; zhi = Domain.zhi;

N1 = options.bimodal.N1;  % average N for type 1 bonds
N2 = options.bimodal.N2;   % average N for type 2 bonds
b  = options.b;

% Cutoff selection 
r1_avg = b*sqrt(N1); 
r2_avg = b*sqrt(N2);


%% Check for advanced options
if (advancedOptions.iadvancedoptions)

    bm = advancedOptions.bimodal;
    % Calculate inital bond pre-stretch (end-to-end distance)
    sig1 = bm.bin_std1_factor*r1_avg;      % standard deviation of bond length distribution
    FWHM1 = bm.bin_width1_factor*sig1;     % full width at half maximum 

    sig2 = bm.bin_std2_factor*r2_avg;
    FWHM2 = bm.bin_width2_factor*sig2;

else

    % use default values
    % Calculate inital bond pre-stretch (end-to-end distance)
    sig1 = 0.4*r1_avg;      % standard deviation of bond length distribution
    FWHM1 = 2.355*sig1;     % full width at half maximum 

    sig2 = 1.5*sig1;
    FWHM2 = 2.355*sig2;

end

% Bond cutoff ranges
r1_lower_cutoff = (r1_avg - FWHM1/2)^2;     %type 1         
r1_upper_cutoff = (r1_avg + FWHM1/2)^2;

r2_lower_cutoff = (r2_avg - FWHM2/2)^2;     %type 2
r2_upper_cutoff = (r2_avg + FWHM2/2)^2;

% some warning checks
if r1_upper_cutoff < 0
    warning("lower cutoff for bonding negative distance: Increase r1_avg or decrease std")
end

if r2_upper_cutoff < 0
    warning("lower cutoff for bonding negative distance: Increase r2_avg or decrease std")
end

if r2_lower_cutoff < r1_upper_cutoff
    warning("overlap in bond length cutoffs between type 1 and type 2 bonds")
end

% Assign method for prob
if strcmp(options.bimodal.distribution_height_mode,'prob')
    P2 = options.bimodal.P;             % desired fraction of type 2 bonds
elseif strcmp(options.bimodal.distribution_height_mode,'fixed')
    N2_number = options.bimodal.N2_number;  % fraction of type

    if N2_number > Max_bond
        warning(" requested number of type 2 bonds exceeds maximum number of bonds")
    end

end

% --------- Unpack Atoms ---------
ids = Atoms(:,1);         % IDs by row
x   = Atoms(:,2);  
y = Atoms(:,3);

% Build ID->row map (robust to id ~= row)
id2row = containers.Map('KeyType','int64','ValueType','int32');
for r = 1:natom
    id2row(int64(ids(r))) = int32(r);
end

% --------- Connect nearby nodes ---------
Bonds      = zeros(Max_bond, 5);
nbond = 0;
no_progress = 0;
ntries = 0;

tic
while (nbond < Max_bond) && (no_progress < stall_limit) %&& (ntries < global_limit)
    
    ntries = ntries + 1;

    % Randomly select distinct node
    r1 = randi(natom);
    if Atoms(r1,5) >= Max_peratom_bond
        no_progress = no_progress + 1;
        continue;
    end

    % return list of possible r2 candidates based on distance cutoffs
    candidate_list = [];
    candidate_types = [];
    for idx = 1:natom

        if r1 == idx, continue; end
        if Atoms(idx,5) >= Max_peratom_bond, continue; end
        dxv = x(idx)-x(r1);
        dyv = y(idx)-y(r1);
        L = dxv*dxv + dyv*dyv;

        if ((L >= r1_lower_cutoff) && (L <= r1_upper_cutoff))
            candidate_list(end+1) = idx; %#ok<AGROW>
            candidate_types(end+1) = 1; %#ok<AGROW>
        end
        if (L >= r2_lower_cutoff) && (L <= r2_upper_cutoff)
            candidate_list(end+1) = idx; %#ok<AGROW>
            candidate_types(end+1) = 2; %#ok<AGROW>
        end

    end

    if isempty(candidate_list)
        %no_progress = no_progress + 1;
        continue;
    end

    % Remove candidates that are already bonds
    for k = 1:Atoms(r1,5)
        bonded_id = Atoms(r1,5 + k);
        bonded_row = find(ids == bonded_id, 1);
        remove_indices = find(candidate_list == bonded_row);
        candidate_list(remove_indices) = [];
        candidate_types(remove_indices) = [];
    end

    if isempty(candidate_list)
        %no_progress = no_progress + 1;
        continue;
    end

    % Select candidate based on probability or fixed number
    if strcmp(options.bimodal.distribution_height_mode,'prob')
        if rand() < P2
            % Randomly select type 2 bond
            type2_indices = find(candidate_types == 2);
            if isempty(type2_indices)
                %no_progress = no_progress + 1;
                continue; % invalid
            end 
            select_idx = type2_indices(randi(length(type2_indices)));
            
        else
            % Randomly select type 1 bond
            type1_indices = find(candidate_types == 1);
            if isempty(type1_indices)
                %no_progress = no_progress + 1;
                continue; % invalid
            end
            select_idx = type1_indices(randi(length(type1_indices)));
            
        end
    elseif strcmp(options.bimodal.distribution_height_mode,'fixed')
        if sum(Bonds(:,5)==2) < N2_number
            % Randomly select type 2 bond
            type2_indices = find(candidate_types == 2);
            if isempty(type2_indices)
                %no_progress = no_progress + 1;
                continue; % invalid
            end 
            select_idx = type2_indices(randi(length(type2_indices)));
            
        else
            % Randomly select type 1 bond
            type1_indices = find(candidate_types == 1);
            if isempty(type1_indices)
                %no_progress = no_progress + 1;
                continue; % invalid
            end
            select_idx = type1_indices(randi(length(type1_indices)));
        end
    end

    r2 = candidate_list(select_idx);
    bond_type = candidate_types(select_idx);
    L  = sqrt((x(r2)-x(r1))^2 + (y(r2)-y(r1))^2);

    % Bond bookkeeping
    nbond = nbond + 1;
    Bonds(nbond,:) = [nbond, ids(r1), ids(r2), sqrt(L), bond_type];
        
    Atoms(r1,5) = Atoms(r1,5) + 1;
    Atoms(r1,5 + Atoms(r1,5)) = ids(r2);

    Atoms(r2,5) = Atoms(r2,5) + 1;
    Atoms(r2,5 + Atoms(r2,5)) = ids(r1);

    no_progress = 0;
end

%stats
type1 = Bonds(:,5) == 1;
type2 = Bonds(:,5) == 2;
Ntotal_kuhn = N1*sum(type1) + N2*sum(type2);
Nkuhn_avg = (N1*sum(type1) + N2*sum(type2))/(sum(type1)+sum(type2));

fprintf('   Placed %d bonds in %4.4f sec \n', nbond, toc);
if strcmp(options.bimodal.distribution_height_mode,'fixed')
    fprintf('   %4.0f type 2 bonds, requested %4.0f \n',sum(type2),N2_number);
else
    fprintf('   %4.2f percent type 2 bonds, requested %4.2f \n',sum(type2)/nbond,P2);
end

if ntries >= global_limit
    warning('Bond creation: hit global try limit (%d).', global_limit);
end
if no_progress >= stall_limit
    warning('Bond creation: local stall after %d attempts.', stall_limit);
end

Bonds = Bonds(1:nbond,:);


% --------- Iterative pruning (row space), no ID compaction ----------
if ~isempty(Bonds)
    changed = true;
    while changed
        % recompute degree from current bonds
        deg = Atoms(:,5);
        %for k=1:size(Bonds,1)
        %    deg(Bonds(k,1)) = deg(Bonds(k,1)) + 1;
        %    deg(Bonds(k,2)) = deg(Bonds(k,2)) + 1;
        %end
        to_del = find(deg <= min_keep);
        if isempty(to_del)
            changed = false;
            break;
        end
        kill = false(size(Bonds,1),1);
        mark = false(natom,1); mark(to_del) = true;
        for k=1:size(Bonds,1)
            if mark(Bonds(k,2)) || mark(Bonds(k,3))
                kill(k) = true;
            end
        end
        if any(kill)
            Bonds = Bonds(~kill,:);
            changed = true;
        else
            changed = false;
        end
    end

    % refresh lengths from coords (robust)
    for k=1:size(Bonds,1)
        r1 = Bonds(k,2); r2 = Bonds(k,3);
        dxv = x(r2)-x(r1); dyv = y(r2)-y(r1);
        Bonds(k,4) = sqrt(dxv*dxv + dyv*dyv);
    end
end


% --------- Build BondsOut in ID space; renumber bond IDs -----------
nb = size(Bonds,1);
BondsOut = zeros(nb,4);
for k=1:nb
    r1 = Bonds(k,2); r2 = Bonds(k,3);
    id1 = ids(r1); id2 = ids(r2);
    L   = Bonds(k,4);
    BondsOut(k,:) = [k, id1, id2, L];
end

% --------- Rebuild Atoms neighbors using **IDs** -------------------
% Clear neighbor metadata
Atoms(:,5) = 0;                             % num_bond
if size(Atoms,2) < 5+Max_peratom_bond
    % extend columns if needed (safeguard)
    Atoms(:, size(Atoms,2)+1 : 5+Max_peratom_bond) = 0;
else
    Atoms(:, 6 : 5+Max_peratom_bond) = 0;   % zero neighbor IDs
end

% Populate neighbor IDs from final bonds
for k=1:nb
    id1 = BondsOut(k,2);  id2 = BondsOut(k,3);
    r1  = id2row(int64(id1));
    r2  = id2row(int64(id2));

    % r1 side
    nb1 = Atoms(r1,5);
    if nb1 < Max_peratom_bond
        Atoms(r1,5) = nb1 + 1;
        Atoms(r1,5 + Atoms(r1,5)) = id2;   % store neighbor **ID**
    end
    % r2 side
    nb2 = Atoms(r2,5);
    if nb2 < Max_peratom_bond
        Atoms(r2,5) = nb2 + 1;
        Atoms(r2,5 + Atoms(r2,5)) = id1;   % store neighbor **ID**
    end
end

% Remove atoms with no bonds
%del = Atoms(:,5) == 0;
%Atoms(del,:) = [];

%ndel = sum(del); natom = natom - ndel;

AtomsOut = Atoms;
fprintf('   Pruned %d atoms, and %d bonds\n',natom-length(AtomsOut),nbond-length(BondsOut));



end