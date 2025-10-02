% -------------------------------------------------------------------------
% Polydisperse mesoscale network generator (2D, no PBC) — MATLAB R2016a
% - Enforces minimum spacing between scattered nodes (rejection sampling)
% - Randomly connects nearby nodes under distance cutoff (with guards)
% - Assigns bimodal bond distribution
% - Iteratively prunes nodes with degree <= 2
% - Exports:
%     * Network.txt  (LAMMPS data: atoms/bonds; atom type=1, bond type=1)
%     * bond.table   (# Chain stats; KEY; N <#bonds>; lines: id i j N b)
% - Per-bond Kuhn segments inferred: N = max(min_N, rounding((L/b)^2))
% -------------------------------------------------------------------------

clc; clear; close all;

%% --------------------------- USER INPUTS --------------------------------
% Domain (2D)
xlo = -10;   xhi =  10;
ylo = -10;   yhi =  10;
zlo =  -1;   zhi =   1;         % thin z just to satisfy LAMMPS format

% Network size caps
Max_atom           = 3200;      % requested # nodes
Max_peratom_bond   = 5;         % degree cap per node (fits neighbor slots)
Max_bond           = 9875;      % cap on bonds

% Bond Khun segments
N1 = 50;
N2 = 250; %[50 350 500]

% Kuhn length and connectivity rule
b = 0.05;                      % constant Kuhn length

% Geometric based N ~ (L/b)^2
r1_avg = b*sqrt(N1);           % average end-to-end bond length for btype 1                    
r2_avg = b*sqrt(N2);           % average end-to-end bond length for btype 2

r2_avg/r1_avg

sig1 = 0.4*r1_avg;             % standard deviation in end-to-end lengths
FWHM1 = 2.355*sig1;            % full width at half maximum 

sig2 = sig1;
FWHM2 = 2.355*sig2;

P2 = 0.1;                    % Probability constraint to reduce number of long bonds

% Cutoff ranges
r1_lower_cutoff = r1_avg - FWHM1/2;              
r1_upper_cutoff = r1_avg + FWHM1/2;

r2_lower_cutoff = r2_avg - FWHM2/2;
r2_upper_cutoff = r2_avg + FWHM2/2;

% some warning checks
if r1_upper_cutoff < 0
    warning("lower cutoff for bonding negative distance: Increase r1_avg or decrease std")
end
% if r1_upper_cutoff > r2_lower_cutoff
%     middistance = (r2_avg - r1_avg) + r1_avg;
%     r1_upper_cutoff = middistance;
%     r2_lower_cutoff = middistance;
% end

% Node scattering (min distance) + guards
min_node_sep               = 0.3;   % HARD minimum spacing between any 2 nodes 0.30
min_node_sep2              = min_node_sep^2;        % compare squared dists
node_scatter_max_tries     = 10 * Max_atom;         % global limit for scatter loop
max_tries_per_node_sample  = 200;                   % per-node rejection limit


connect_cutoff     = 0.7 * b * 40;          % distance to allow bond creation

% Bond creation guards
bond_global_try_limit          = 200 * Max_bond; % absolute ceiling
max_attempts_without_progress  = 10  * Max_atom; % local stall guard

% Pruning rule
min_degree_keep    = 2;        % delete nodes with degree <= 2

% Kuhn inference
kuhn_rounding      = 'round';  % 'round' | 'ceil' | 'floor'
min_N              = 1;        % minimum segments per bond

% Kuhn inference mode
% 'geom'  -> infer from L via N = (L/b)^2 (then rounded)
% 'range' -> map lengths monotonically into [N_target_min, N_target_max]
N_assignment_mode = 'range';      % 'geom' or 'range'
N_range_method    = 'rank';    % 'linear' or 'rank' (only used if mode='range')
N_target_min      = 10;          % min N when mode='range'
N_target_max      = 80;          % max N when mode='range'

% RNG seed (fix for reproducibility; comment to randomize each run)
rng(1);

% Output files
lammps_data_file   = 'PronyNetwork.txt';
lammps_visual_file = 'PronyVisual.dat';
bond_table_file    = 'bond.table';

%% ---------------------- SCATTER NODES (WITH MIN DIST) -------------------
% Atoms layout:
% [ ID | X | Y | Z | num_bond | nbr1 | nbr2 | nbr3 | nbr4 | spare ]
Atoms  = zeros(Max_atom, 10);
N_atom = 0;

global_scatter_tries = 0;

while (N_atom < Max_atom) && (global_scatter_tries < node_scatter_max_tries)
    global_scatter_tries = global_scatter_tries + 1;

    % Try to draw a point that respects min_node_sep from all existing nodes
    accepted = false;
    per_node_tries = 0;

    while (~accepted) && (per_node_tries < max_tries_per_node_sample)
        per_node_tries = per_node_tries + 1;

        xi = xlo + (xhi - xlo) * rand;
        yi = ylo + (yhi - ylo) * rand;
        zi = 0; % planar

        if N_atom == 0
            accepted = true;  % first node always accepted
        else
            % Check squared distances to all currently placed nodes
            ok = true;
            for k = 1:N_atom
                dx = xi - Atoms(k,2);
                dy = yi - Atoms(k,3);
                if (dx*dx + dy*dy) < min_node_sep2
                    ok = false;
                    break;
                end
            end
            accepted = ok;
        end
    end

    if ~accepted
        % Could not find a valid spot for this node within per-node tries;
        % move on (global loop will stop when global_scatter_tries hits cap).
        continue;
    end

    % Place node
    N_atom = N_atom + 1
    Atoms(N_atom,1) = N_atom;
    Atoms(N_atom,2) = xi;
    Atoms(N_atom,3) = yi;
    Atoms(N_atom,4) = zi;
end

if N_atom < Max_atom
    warning('Requested %d atoms, placed %d (scatter tries=%d / %d).', ...
            Max_atom, N_atom, global_scatter_tries, node_scatter_max_tries);
    if N_atom == 0
        error('No atoms placed—aborting.');
    end
    Atoms = Atoms(1:N_atom, :);
end

figure; scatter(Atoms(1:N_atom,2), Atoms(1:N_atom,3), 6, 'k', 'filled');
axis equal tight; title(sprintf('Scattered nodes (min sep = %.3g)', min_node_sep));

%% ---------------------- CONNECT NEARBY NODES ----------------------------
% Bonds layout: [ bondID | i1 | i2 | L | type]
Bonds      = zeros(Max_bond, 5);
Total_bond = 0;

attempts_without_progress = 0;
bond_global_try_count     = 0;

while (Total_bond < Max_bond) && ...
      (attempts_without_progress < max_attempts_without_progress) && ...
      (bond_global_try_count < bond_global_try_limit)

    bond_global_try_count = bond_global_try_count + 1;

    % Pick a random unsaturated node
    i1 = randi(N_atom);
    if (Atoms(i1,1) == 0) || (Atoms(i1,5) >= Max_peratom_bond)
        attempts_without_progress = attempts_without_progress + 1;
        continue;
    end

    x1 = Atoms(i1,2); y1 = Atoms(i1,3);
    
    % Build neighbor candidate list (within cutoff, unsaturated, not connected)
    Num_neigh  = 0;
    Neigh_list = zeros(1, 5); % grows as needed
    Type_list = zeros(1,5);

    nb_i1 = Atoms(i1,5);

    for j = 1:N_atom
        if j == i1, continue; end
        if (Atoms(j,1) == 0) || (Atoms(j,5) >= Max_peratom_bond), continue; end

        dx = Atoms(j,2) - x1;
        dy = Atoms(j,3) - y1;
        Dist = sqrt(dx*dx + dy*dy);

        %% See if it meets requirements to be type 1
        if r1_lower_cutoff < Dist && Dist < r1_upper_cutoff
            % Already connected? (scan neighbor slots of i1)
            is_connected = 0;
            if nb_i1 > 0
                for ss = 1:nb_i1
                    if Atoms(i1,5 + ss) == j
                        is_connected = 1;
                        break;
                    end
                end
            end
            if ~is_connected
                Num_neigh = Num_neigh + 1;
                if Num_neigh > numel(Neigh_list)
                    tmp = zeros(1, numel(Neigh_list) + 50);
                    tmp(1:numel(Neigh_list)) = Neigh_list;
                    Neigh_list = tmp; % grow array (R2016a-safe)

                    tmp_ty = zeros(1, numel(Type_list) + 5);
                    tmp_ty(1:numel(Type_list)) = Type_list;
                    Type_list = tmp_ty;
                end
                Neigh_list(Num_neigh) = j;
                Type_list(Num_neigh) = 1;
            end
        end

        %%% See if it meets requirements to be type 2
        if r2_lower_cutoff < Dist && Dist < r2_upper_cutoff
            % Already connected? (scan neighbor slots of i1)
            is_connected = 0;
            if nb_i1 > 0
                for ss = 1:nb_i1
                    if Atoms(i1,5 + ss) == j
                        is_connected = 1;
                        break;
                    end
                end
            end
            if ~is_connected
                Num_neigh = Num_neigh + 1;
                if Num_neigh > numel(Neigh_list)
                    tmp = zeros(1, numel(Neigh_list) + 5);
                    tmp(1:numel(Neigh_list)) = Neigh_list;
                    Neigh_list = tmp; % grow array (R2016a-safe)

                    tmp_ty = zeros(1, numel(Type_list) + 5);
                    tmp_ty(1:numel(Type_list)) = Type_list;
                    Type_list = tmp_ty;
                end
                Neigh_list(Num_neigh) = j;
                Type_list(Num_neigh) = 2;
            end
        end
    end

    if Num_neigh == 0
        attempts_without_progress = attempts_without_progress + 1;
        continue;
    end
    
    % Draw a probability to determine what bond type to attempt
    iP2 = rand() < P2;

    if iP2
        idelete = Type_list == 1; % delete type1's from list
        Neigh_list(idelete) = [];
        Type_list(idelete) = [];
        Num_neigh = Num_neigh - sum(idelete);
    else
        idelete = Type_list == 2; % delete type2's from list
        Neigh_list(idelete) = [];
        Type_list(idelete) = [];
        Num_neigh = Num_neigh - sum(idelete);
    end

    if Num_neigh == 0
        attempts_without_progress = attempts_without_progress + 1;
        continue;
    end

    % Pick a neighbor and form a bond
    id_random = randi(Num_neigh);
    i2 = Neigh_list( id_random );
    type = Type_list( id_random) ;
    if (Atoms(i2,5) >= Max_peratom_bond)
        attempts_without_progress = attempts_without_progress + 1;
        continue;
    end

    % neighbor bookkeeping (add i2 to i1, and i1 to i2)
    Atoms(i1,5) = Atoms(i1,5) + 1;
    Atoms(i1,5 + Atoms(i1,5)) = i2;

    Atoms(i2,5) = Atoms(i2,5) + 1;
    Atoms(i2,5 + Atoms(i2,5)) = i1;

    % record bond
    Total_bond = Total_bond + 1;

    x2 = Atoms(i2,2); y2 = Atoms(i2,3);
    L  = sqrt((x2 - x1)^2 + (y2 - y1)^2);

    Bonds(Total_bond,1) = Total_bond;
    Bonds(Total_bond,2) = i1;
    Bonds(Total_bond,3) = i2;
    Bonds(Total_bond,4) = L;
    Bonds(Total_bond,5) = type;

    attempts_without_progress = 0; % reset on success
end

if bond_global_try_count >= bond_global_try_limit
    warning('Stopped bond creation: hit global try limit (%d).', bond_global_try_limit);
end
if attempts_without_progress >= max_attempts_without_progress
    warning('Stopped bond creation: local stall (%d attempts w/o progress).', ...
            max_attempts_without_progress);
end

% Plot bonds (pre-prune)
figure; hold on;
scatter(Atoms(1:N_atom,2), Atoms(1:N_atom,3), 6, 'k', 'filled');
for k = 1:Total_bond
    if Bonds(k,1) == 0, continue; end
    i1 = Bonds(k,2); i2 = Bonds(k,3);
    if Bonds(k,5) == 1
        plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'k-');
    else
        plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'r-');
    end
    
end
axis equal tight; title('Initial bonds (pre-prune)');

%% ---------------------- PRUNE LOW-DEGREE NODES --------------------------
% Iteratively delete nodes with degree <= 2 and their incident bonds
changed = true;
while changed
    changed = false;

    to_delete = (Atoms(:,1) ~= 0) & (Atoms(:,5) <= (min_degree_keep - 1));
    if any(to_delete)
        changed = true;
        ids_to_delete = Atoms(to_delete,1);

        % Remove their references from neighbors and decrement counts
        for idx = 1:numel(ids_to_delete)
            Must_Del = ids_to_delete(idx);

            % For every atom, if neighbor == Must_Del, remove it
            for ii = 1:N_atom
                if Atoms(ii,1) == 0, continue; end
                nb = Atoms(ii,5);
                if nb <= 0, continue; end

                slot = 1;
                while slot <= nb
                    if Atoms(ii,5 + slot) == Must_Del
                        % Move last neighbor into this slot
                        if slot ~= nb
                            Atoms(ii,5 + slot) = Atoms(ii,5 + nb);
                        end
                        Atoms(ii,5 + nb) = 0;
                        Atoms(ii,5) = Atoms(ii,5) - 1;
                        nb = nb - 1;
                        % do not increment slot (new value now at this slot)
                    else
                        slot = slot + 1;
                    end
                end
            end

            % Remove bonds that connect to Must_Del
            for bb = 1:Total_bond
                if Bonds(bb,1) == 0, continue; end
                if (Bonds(bb,2) == Must_Del) || (Bonds(bb,3) == Must_Del)
                    Bonds(bb,:) = 0;
                end
            end

            % Zero-out atom row
            Atoms(Must_Del,:) = 0;
        end
    end
end

% Compact Atoms -> All_Atoms (preserve original IDs)
All_Atoms = zeros(sum(Atoms(:,1)~=0), size(Atoms,2));
Atom_count = 0;
for i = 1:N_atom
    if Atoms(i,1) ~= 0
        Atom_count = Atom_count + 1;
        All_Atoms(Atom_count,:) = Atoms(i,:);
    end
end

% Compact Bonds -> All_Bonds; drop any referencing deleted atoms
Bond_count = 0;
All_Bonds = zeros(0,5);
if Total_bond > 0
    tmpB = zeros(Total_bond, 5);
    for n = 1:Total_bond
        if Bonds(n,1) == 0, continue; end
        i1 = Bonds(n,2); i2 = Bonds(n,3);
        if (i1==0) || (i2==0) || (Atoms(i1,1)==0) || (Atoms(i2,1)==0)
            continue;
        end
        Bond_count = Bond_count + 1;
        tmpB(Bond_count,:) = Bonds(n,:);
    end
    All_Bonds = tmpB(1:Bond_count, :);
end

% Recompute bond lengths from surviving geometry
for k = 1:Bond_count
    i1 = All_Bonds(k,2);
    i2 = All_Bonds(k,3);
    x1 = Atoms(i1,2); y1 = Atoms(i1,3);
    x2 = Atoms(i2,2); y2 = Atoms(i2,3);
    All_Bonds(k,4) = sqrt((x2-x1)^2 + (y2-y1)^2);
end

% Visual check after pruning
figure; hold on;
scatter(All_Atoms(:,2), All_Atoms(:,3), 6, 'k', 'filled');
for k = 1:Bond_count
    i1 = All_Bonds(k,2); i2 = All_Bonds(k,3);
    if All_Bonds(k,5) == 1
        plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'k-');
    else
        plot([Atoms(i1,2) Atoms(i2,2)], [Atoms(i1,3) Atoms(i2,3)], 'r-');
    end
end
axis equal tight; title('Final network (post-prune)');

% Quick histograms
if Bond_count > 0
    figure; hist(All_Bonds(:,4), 80);
    xlabel('Bond length L'); ylabel('Count'); title('Length distribution (final)');
end

% Bimodal histograms
type1 = All_Bonds(:,5) == 1;
figure; hold on; 
histogram(All_Bonds(type1,4),60,'FaceColor',[0.2 0.2 0.2],'FaceAlpha',1,'LineWidth',0.0005);
histogram(All_Bonds(~type1,4),60,'FaceColor',[1 0 0],'FaceAlpha',1,'LineWidth',0.0005);
axis tight
axis([0.3 1.3 0 350]); 
xlabel('Bond length L'); ylabel('Count'); title('Length distribution (final)')
set(gca,'FontSize',16,'LineWidth',2)

fprintf("%4.2f percent type 1 bonds | Desired percent: %4.2f \n",100*sum(type1)/length(type1),100*(1-P2))

%% ---------------------- WRITE LAMMPS DATA: Network.txt ------------------
fid = fopen(lammps_data_file,'w');
if fid < 0, error('Could not open %s for writing.', lammps_data_file); end

fprintf(fid, '\n\n');
fprintf(fid, '%d atoms\n', Atom_count);
fprintf(fid, '%d bonds\n', Bond_count);
fprintf(fid, '%d atom types\n', 1);
fprintf(fid, '%d bond types\n', 2);
fprintf(fid, '%f %f xlo xhi\n', xlo, xhi);
fprintf(fid, '%f %f ylo yhi\n', ylo, yhi);
fprintf(fid, '%f %f zlo zhi\n', zlo, zhi);
fprintf(fid, '\n');

fprintf(fid, 'Atoms #bpm/sphere\n\n');
for i = 1:Atom_count
    % Format: atomID  molID  atomType  x y z
    fprintf(fid, '%d 1 1 1 1 %f %f %f\n', All_Atoms(i,1), All_Atoms(i,2), All_Atoms(i,3), All_Atoms(i,4));
end

fprintf(fid, '\nBonds\n\n');
for i = 1:Bond_count
    % Format: bondID  bondType  atom1  atom2
    fprintf(fid, '%d %d %d %d\n', All_Bonds(i,1), All_Bonds(i,5), All_Bonds(i,2), All_Bonds(i,3));
end

fclose(fid);
fprintf('Wrote %s with %d atoms and %d bonds.\n', lammps_data_file, Atom_count, Bond_count);

%% ------------ WRITE LAMMPS DATA: (visualization) Network.txt ------------
fid = fopen(lammps_visual_file,'w');
if fid < 0, error('Could not open %s for writing.', lammps_visual_file); end

fprintf(fid, '\n\n');
fprintf(fid, '%d atoms\n', Atom_count);
fprintf(fid, '%d bonds\n', Bond_count);
fprintf(fid, '%d atom types\n', 1);
fprintf(fid, '%d bond types\n', 2);
fprintf(fid, '%f %f xlo xhi\n', xlo, xhi);
fprintf(fid, '%f %f ylo yhi\n', ylo, yhi);
fprintf(fid, '%f %f zlo zhi\n', zlo, zhi);
fprintf(fid, '\n');

% Atom type: hybrid bond sphere
fprintf(fid, 'Atoms #hybrid bond sphere\n\n');
for i = 1:Atom_count
    % Format: atom-ID atom-type x y z molecule-ID diameter density
    fprintf(fid, '%d 1 %f %f %f 1 1 1 \n', All_Atoms(i,1), All_Atoms(i,2), All_Atoms(i,3), All_Atoms(i,4));
end

fprintf(fid, '\nBonds\n\n');
for i = 1:Bond_count
    % Format: bondID  bondType  atom1  atom2
    fprintf(fid, '%d %d %d %d\n', All_Bonds(i,1), All_Bonds(i,5), All_Bonds(i,2), All_Bonds(i,3));
end

fclose(fid);
fprintf('Wrote %s with %d atoms and %d bonds.\n', lammps_visual_file, Atom_count, Bond_count);

%% ---------------------- WRITE bond.table (polydisperse) -----------------
%% ---------------------- WRITE bond.table (polydisperse) -----------------
% if Bond_count > 0
%     Lvec = All_Bonds(:,4);
% 
%     % Compute Nvec by selected mode
%     if strcmpi(N_assignment_mode,'geom')
%         % Geometry-based: N ~ (L/b)^2
%         switch lower(kuhn_rounding)
%             case 'ceil'
%                 Nvec = ceil( (Lvec./b).^2 );
%             case 'floor'
%                 Nvec = floor( (Lvec./b).^2 );
%             otherwise % 'round'
%                 Nvec = round( (Lvec./b).^2 );
%         end
%         Nvec = max(Nvec, min_N);
% 
%     else
%         % Range-mapped: shortest bond -> N_target_min, longest -> N_target_max
%         % Sanity & ordering
%         Nlo = min(N_target_min, N_target_max);
%         Nhi = max(N_target_min, N_target_max);
% 
%         if strcmpi(N_range_method,'rank')
%             % Rank-based monotone assignment with uniform coverage across bonds
%             [Ls, idx] = sort(Lvec, 'ascend');
%             % Create evenly spaced targets across [Nlo, Nhi]
%             if Bond_count == 1
%                 Ntargets = (Nlo + Nhi)/2;
%             else
%                 Ntargets = linspace(Nlo, Nhi, Bond_count);
%             end
%             Ntmp = zeros(Bond_count,1);
%             for k = 1:Bond_count
%                 Ntmp(idx(k)) = Ntargets(k);
%             end
%             Nvec = round(Ntmp);
%         else
%             % Linear scaling by actual length range
%             Lmin = min(Lvec);  Lmax = max(Lvec);
%             if Lmax == Lmin
%                 Nvec = round( (Nlo + Nhi)/2 * ones(Bond_count,1) );
%             else
%                 t = (Lvec - Lmin) ./ (Lmax - Lmin);   % in [0,1]
%                 Nvec = round( Nlo + t .* (Nhi - Nlo) );
%             end
%         end
% 
%         % Clamp and enforce at least min_N
%         Nvec = max(Nvec, min_N);
%         Nvec = min(Nvec, max(N_target_min, N_target_max));
%     end
% 
%     % Quick sanity plot of N distribution
%     figure; hist(Nvec, max(10, min(80, ceil(sqrt(numel(Nvec))))));
%     xlabel('Assigned N per bond'); ylabel('Count');
%     title(sprintf('N distribution (%s%s)', ...
%         N_assignment_mode, ...
%         strcmpi(N_assignment_mode,'range'), N_range_method));
% 
%     % Write the table
%     fidBT = fopen(bond_table_file,'w');
%     if fidBT < 0, error('Could not open %s for writing.', bond_table_file); end
% 
%     fprintf(fidBT, '# Chain stats\n\n');
%     fprintf(fidBT, 'KEY\n');
%     fprintf(fidBT, 'N %d\n\n', Bond_count);
% 
%     % Lines: bond-id   i1-atom   i2-atom   N    b
%     for k = 1:Bond_count
%         fprintf(fidBT, '%d %d %d %d %.8g\n', ...
%             All_Bonds(k,1), All_Bonds(k,2), All_Bonds(k,3), Nvec(k), b);
%     end
% 
%     fclose(fidBT);
%     fprintf('Wrote %s with %d entries. b = %.6g, mode=%s\n', ...
%         bond_table_file, Bond_count, b, N_assignment_mode);
% else
%     warning('No bonds to write into %s.', bond_table_file);
% end

