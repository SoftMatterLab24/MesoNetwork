% -------------------------------------------------------------------------
% Creates a single edge notch in an equilibriated network
% - Finds bonds that intersect crack vector and removes from bondlist
% - Crack length normalized w.r.t. network mesh size (user-specified)
% - 
% - 
% - Exports:
%    * NetworkCrack.txt with updated bond list
% -------------------------------------------------------------------------

clc; clear; close all;

%% --------------------------- USER INPUTS --------------------------------

xi = 1.6; %mesh size (LJ lengthscale)
crack_mode = 'notch'; % 'sharp' or 'notch'

%% Crack dimensions
c       = 2000*xi;     % notch length
t       = 75*xi;       % notch thickness
alpha   = 90;         % taper angle (deg)

%!!! MUST BE UPDATED !!!
xmin = -2540; %leftmost xboundary of the network (look at .dump)

%location to lammps outputs (atoms_equilib.dump,bonds_equilib.dump)
% <<<<<<< HEAD
loc = 'G:\LAMMPS_data\Notch_length_study\EX_MD_smp1';
% =======
% loc = 'C:\Users\zwhit\Downloads\polydisperse_net_generator\Runs\bimodal\Unnotched\002';
% >>>>>>> 751163b0431be0507b5c99b0adcda654be0fa46b

%dump names
atom_name = 'atoms_equilib.dump';
bond_name = 'bonds_equilib.dump';
table_name ='bond.table';

ivisual = 1; %create visual dump

%new write names
lammps_data_file   = 'PolyNetwork_MD_SMP0001_N0001_crack.dat';
lammps_table_file  = 'bondC_MD_SMP0001_N0001.table';
lammps_visual_file = 'PronyNetworkCrack_VISUAL2_1000.dat';
%% --------------------------- Create Crack -------------------------------

% filenames for dumps
afile = strcat(loc,'\',atom_name);
bfile = strcat(loc,'\',bond_name);
tfile = strcat(loc,'\',table_name);

%parse data from equilibrated neworks
[x,y,id,ty,mol,natoms,bAtom1,bAtom2,bType,bLength,bForce,xlims,ylims,zlims] = parse_equidump_full_fun(afile,bfile);

%parse data from bond.table
[id_t, iatom_t, jatom_t, N, b] = parse_bond_table(tfile);

%%% Extract data from the last time step

%limits
Xlimits = xlims{end};
Ylimits = ylims{end};
Zlimits = zlims{end};

dX = Xlimits(2) - Xlimits(1);
dY = Ylimits(2) - Ylimits(1);

%coords
X = dX*x{end}-Xlimits(2);
Y = dY*y{end}-Ylimits(2);

%types
TYPE = ty{end};
bondType = bType{end};

%typology
iatom = bAtom1{end};
jatom = bAtom2{end};

%IDs
ID = id{end};
MOL = mol{end};

%get bond xcentroid
%for ii = 1:length(bondType)
%    xlmin = [X(iatom(ii)) X(jatom(ii))];
%    xcen(ii) = (X(iatom(ii)) + X(jatom(ii)))/2;
%    ycen(ii) = (Y(iatom(ii)) + Y(jatom(ii)))/2;
%end

xmax = max(X);
ymin = min(Y);
ymax = max(Y);
dely = (ymax-ymin)/2;

%%% Clean bond table if it somehow has extra entries (i.e. bonds broke during sim) %%%
if length(iatom_t) ~= length(iatom)
    
    ikeep = logical(zeros([length(iatom_t) 1]));
    % loop through bonds in table
    for i = 1:length(id_t)
        itab = iatom_t(i); jtab = jatom_t(i);

        for j = 1:length(iatom)
           if (itab == iatom(j) && jtab == jatom(j)) || (jtab == iatom(j) && itab == jatom(j)) 
                ikeep(i) = true;
           end
        end
    end

    %clean table
    id_t(~ikeep)    = [];
    iatom_t(~ikeep) = [];
    jatom_t(~ikeep) = [];
    N(~ikeep)       = [];
    b(~ikeep)       = [];
end

%%% remove bonds or atoms based on crack type
if strcmpi(crack_mode,'sharp')
    % Sharp crack
    sortedIDs = NaN(length(ID),1);
    for ii = 1:length(sortedIDs)
        sortedIDs(ID(ii)) = ii;
    end
    
    i_index = sortedIDs(iatom(:));
    j_index = sortedIDs(jatom(:));
        
    %list of bond vectors
    XY1 = [X(i_index) Y(i_index) X(j_index) Y(j_index)];

    %crack vector
    XY2 = [Xlimits(1) dely+ymin xmin+c dely+ymin];

    %find where crack intersects bonds
    out = lineSegmentIntersect(XY1,XY2);
    idx = out.intAdjacencyMatrix;

    %update
    Atom_count = length(ID);
    Bond_count = sum(~idx);

    ikeep = logical(zeros([length(idx) 1]));
    % loop through bonds
    for i = 1:length(idx)
        ia = iatom(i); ja = jatom(i);

        %find the index in the bond table
        for j = 1:length(id_t)
            if ((iatom_t(j) == ia && jatom_t(j) == ja) || (iatom_t(j) == ja && jatom_t(j) == ia)) && ~idx(i)
                ikeep(i) = true;
            end
        end
    end

    % remove marked bonds
    id_t(~ikeep)    = [];
    iatom_t(~ikeep) = [];
    jatom_t(~ikeep) = [];
    N(~ikeep)       = [];
    b(~ikeep)       = [];

    % renumber bond IDs
    for ii = 1:length(id_t)
        id_t(ii) = ii;
    end

    iatomn = iatom;
    jatomn = jatom;

    iatom_t_new = iatom_t;
    jatom_t_new = jatom_t;
    
elseif strcmpi(crack_mode,'notch')

    sortedIDs = NaN(length(ID),1);
    for ii = 1:length(sortedIDs)
        sortedIDs(ID(ii)) = ii;
    end
    
    i_index = sortedIDs(iatom(:));
    j_index = sortedIDs(jatom(:));
    
    idx_a   = zeros([1 length(ID)]);
    idx     = zeros([1 length(i_index)]);

    % compute bond centroids
    %xcen = (X(i_index) + X(j_index))/2;
    %ycen = (Y(i_index) + Y(j_index))/2;

    % Find atoms in notch region
    tip_dx = (t/2)*(tand(alpha/2));

    notch_xmin = 1.1*xmin;
    notch_xmax = notch_xmin + c - tip_dx;

    notch_ymin = -t/2;
    notch_ymax = t/2;

    % rectangular region first
    inotchtmp_atoms = find(X >= notch_xmin & X <= notch_xmax & Y >= notch_ymin & Y <= notch_ymax); % atoms in rectangular notch region
    %inotchtmp_bonds = find(xcen >= notch_xmin & xcen <= notch_xmax & ycen >= notch_ymin & ycen <= notch_ymax); % bonds in rectangular notch region
    
    % triangular tip region second
    tip_xmin = notch_xmax;
    tip_xmax = notch_xmax + tip_dx;
    tip_ymax = t/2;
    tip_ymin = -t/2;

    % for atoms
    for ii = 1:length(X)
        if X(ii) >= tip_xmin && X(ii) <= tip_xmax
            %calculate ybound at this x location
            m = (0-tip_ymax)/(tip_xmax - tip_xmin); intercept = - m*tip_xmax;
            ybound = m*X(ii) + intercept;

            if Y(ii) <= ybound && Y(ii) >= -ybound
                inotchtmp_atoms = [inotchtmp_atoms; ii];
            end
        end
    end

    idx_a(inotchtmp_atoms) = 1;
    Atom_count = sum(~idx_a);

    %%% Remap the atom ids %%%

    % renumber atom IDs after removing atoms in notch
    sortedIDs_new = NaN(length(ID),1);
    counter = 0;
    for ii = 1:length(ID)
        if ~idx_a(ii)
            counter = counter + 1;
            sortedIDs_new(ii) = counter;
        end
    end
    
    %remap the iatom/jatom indexes with new atom IDs
    iatomn = sortedIDs_new(i_index);
    jatomn = sortedIDs_new(j_index);
    
    % Remap bond table indexes
    iatom_t_new = NaN([1 length(iatom_t)]);
    jatom_t_new = iatom_t_new;

    ikeep = zeros([length(id_t) 1]);
    % loop through bond table
    for i = 1:length(id_t)
        ia = iatom_t(i); ja = jatom_t(i);

       %find the index in bond dump
        for j = 1:length(idx)
            if ((iatom(j) == ia && jatom(j) == ja) || (iatom(j) == ja && jatom(j) == ia)) % && ~idx(i)
                ikeep(i) = j;
                iatom_t_new(i) = iatomn(j); %get remapped indexes
                jatom_t_new(i) = jatomn(j);
            end
        end
    end
    
    %%% Now that things are remapped -> proceed to delete %%%
    
    % mark bonds for deletion if they are missing an atom (i.e. an atom 
    % index is NaN)
    idx = zeros([length(iatom) 1]);
    for ii = 1:length(iatom)
        if isnan(iatomn(ii)) || isnan(jatomn(ii))
            idx(ii) = true; % Mark bond for deletion
        end
    end
    
    Bond_count = sum(~idx); %update bond count
    
    % mark bond in bond table for deletion
    ikeep2 = logical(ones([length(iatom_t_new) 1]));
    for i = 1:length(iatom_t_new)
        if isnan(iatom_t_new(i)) || isnan(jatom_t_new(i))
            ikeep2(i) = false;
        end
    end
    
    % remove marked bonds
    id_t(~ikeep2)        = [];
    iatom_t_new(~ikeep2) = [];
    jatom_t_new(~ikeep2) = [];
    N(~ikeep2)           = [];
    b(~ikeep2)           = [];

    % renumber bond IDs (for table)
    for ii = 1:length(id_t)
        id_t(ii) = ii;
    end

end

%%% Sanity check to see if bond vector crosses network 
% visual check
figure; hold on
scatter(X,Y,'r.')
scatter(X(inotchtmp_atoms),Y(inotchtmp_atoms),'k.')
% Plot the deleted bonds (faster than all)
for ii = 1:length(iatom)
    if idx(ii)
        ia =  i_index(ii); ja =  j_index(ii);
        plot([X(ia) X(ja)],[Y(ia) Y(ja)],'-k')
    end

end

%% ---------------------- WRITE LAMMPS DATA: Network.txt ------------------

fid = fopen(lammps_data_file,'w');
if fid < 0, error('Could not open %s for writing.', lammps_data_file); end

fprintf(fid, '\n\n');
fprintf(fid, '%d atoms\n', Atom_count);
fprintf(fid, '%d bonds\n', Bond_count);
fprintf(fid, '%d atom types\n', 1);
fprintf(fid, '%d bond types\n', 1);
fprintf(fid, '%f %f xlo xhi\n', 1.01*Xlimits(1), 1.01*Xlimits(2));
fprintf(fid, '%f %f ylo yhi\n', 1.01*Ylimits(1), 1.01*Ylimits(2));
fprintf(fid, '%f %f zlo zhi\n', Zlimits(1), Zlimits(2));
fprintf(fid, '\n');

fprintf(fid, 'Atoms \n\n');
counter = 0;
for i = 1:length(idx_a)
    if ~idx_a(i)
    counter = counter + 1
    % Format: atomID  molID  atomType  x y z
    fprintf(fid, '%d %d %d 1 1 %f %f %f\n', counter, MOL(i), TYPE(i), X(i),Y(i),0);
    end
end

fprintf(fid, '\nBonds\n\n');
counter = 0;
for i = 1:length(idx)
    if ~idx(i)
        counter = counter + 1
        % Format: bondID  bondType  atom1  atom2
        fprintf(fid, '%d %d %d %d\n',counter, bondType(i), iatomn(i),jatomn(i));
    end
end

fclose(fid);
fprintf('Wrote %s with %d atoms and %d bonds.\n', lammps_data_file, Atom_count, Bond_count);

%% ------------ REWRITE Bond.table ------------
fidt = fopen(lammps_table_file,'w');
fprintf(fidt,'Chain statistics\n\n');
fprintf(fidt,'KEY\n');
fprintf(fidt,'N %d\n\n', length(id_t));

for i = 1:length(id_t)
    i
    fprintf(fidt, '%d %d %d %f %f\n', id_t(i), iatom_t_new(i), jatom_t_new(i), N(i), b(i));
end
fclose(fidt);

%% ------------ WRITE LAMMPS DATA: (visualization) Network.txt ------------
if ivisual

fid = fopen(lammps_visual_file,'w');
if fid < 0, error('Could not open %s for writing.', lammps_visual_file); end

fprintf(fid, '\n\n');
fprintf(fid, '%d atoms\n', Atom_count);
fprintf(fid, '%d bonds\n', Bond_count);
fprintf(fid, '%d atom types\n', 1);
fprintf(fid, '%d bond types\n', 1);
fprintf(fid, '%f %f xlo xhi\n', Xlimits(1), Xlimits(2));
fprintf(fid, '%f %f ylo yhi\n', Ylimits(1), Ylimits(2));
fprintf(fid, '%f %f zlo zhi\n', Zlimits(1), Zlimits(2));
fprintf(fid, '\n');

% Atom type: hybrid bond sphere
fprintf(fid, 'Atoms \n\n');
counter = 0;
for i = 1:length(idx_a)
    if ~idx_a(i)
    counter = counter + 1
    % Format: atom-ID atom-type x y z molecule-ID diameter density
    fprintf(fid, '%d 1 %f %f %f 1 1 1 \n', counter, X(i), Y(i), 0);
    end
end

fprintf(fid, '\nBonds\n\n');
counter = 0;
for i = 1:length(idx)
    if ~idx(i)
        counter = counter + 1
        % Format: bondID  bondType  atom1  atom2
        fprintf(fid, '%d %d %d %d\n', counter, bondType(i), iatomn(i), jatomn(i));
    end
end

fclose(fid);
fprintf('Wrote %s with %d atoms and %d bonds.\n', lammps_visual_file, Atom_count, Bond_count);

end