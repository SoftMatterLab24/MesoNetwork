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

xi = 1.6; %mesh size (lengthscale = set to Kuhn length)

%inital notch size
c = 200*xi;

%!!! MUST BE UPDATED !!!
xmin = -480; %leftmost xboundary of the network (look at .dump)

%location to lammps outputs (atoms_equilib.dump,bonds_equilib.dump)
loc = 'C:\Users\zwhit\Downloads\polydisperse_net_generator';

%dump names
atom_name = 'atoms_equilib.dump';
bond_name = 'bonds_equilib.dump';
table_name ='bond.table';

ivisual = 0; %create visual dump

%new write names
lammps_data_file   = 'PolyNetworkCrack_nano_200x200.dat';
lammps_table_file  = 'bondC.table';
%lammps_visual_file = 'PronyNetworkCrack_VISUAL2_1000.dat';
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
for ii = 1:length(bondType)
    xlmin = [X(iatom(ii)) X(jatom(ii))];
    xcen(ii) = (X(iatom(ii)) + X(jatom(ii)))/2;
    %ycen(ii) = (Y(iatom(ii)) + Y(jatom(ii)))/2;
end

xmax = max(X);
ymin = min(Y);
ymax = max(Y);
dely = (ymax-ymin)/2;

%%% create vector for bonds
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

%find centroids
% for ii = 1:length(sortedIDs)
%     cx = (X(i_index) + X(j_index))/2;
% end
% 
% iedge = cx > 0.137;

%update
Atom_count = length(ID);
Bond_count = sum(~idx);

%% Update bond table
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

%% Sanity check to see if bond vector crosses network 
    figure; hold on
    for mm = 1:length(XY1)
        plot([XY1(mm,1) XY1(mm,3)],[XY1(mm,2) XY1(mm,4)])
    end

    plot([XY2(1,1) XY2(1,3)],[XY2(1,2) XY2(1,4)],'k','LineWidth',4)

%% ---------------------- WRITE LAMMPS DATA: Network.txt ------------------
counter = 0;

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
for i = 1:Atom_count
    % Format: atomID  molID  atomType  x y z
    fprintf(fid, '%d %d %d 1 1 %f %f %f\n', ID(i), MOL(i), TYPE(i), X(i),Y(i),0);
end

fprintf(fid, '\nBonds\n\n');
for i = 1:length(idx)
    if ~idx(i)
        counter = counter + 1;
        % Format: bondID  bondType  atom1  atom2
        fprintf(fid, '%d %d %d %d\n',counter, bondType(i), iatom(i),jatom(i));
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
    fprintf(fidt, '%d %d %d %f %f\n', id_t(i), iatom_t(i), jatom_t(i), N(i), b(i));
end
fclose(fidt);

%% ------------ WRITE LAMMPS DATA: (visualization) Network.txt ------------
if ivisual
counter = 0;

fid = fopen(lammps_visual_file,'w');
if fid < 0, error('Could not open %s for writing.', lammps_visual_file); end

fprintf(fid, '\n\n');
fprintf(fid, '%d atoms\n', Atom_count);
fprintf(fid, '%d bonds\n', Bond_count);
fprintf(fid, '%d atom types\n', 1);
fprintf(fid, '%d bond types\n', 2);
fprintf(fid, '%f %f xlo xhi\n', Xlimits(1), Xlimits(2));
fprintf(fid, '%f %f ylo yhi\n', Ylimits(1), Ylimits(2));
fprintf(fid, '%f %f zlo zhi\n', Zlimits(1), Zlimits(2));
fprintf(fid, '\n');

% Atom type: hybrid bond sphere
fprintf(fid, 'Atoms \n\n');
for i = 1:Atom_count
    % Format: atom-ID atom-type x y z molecule-ID diameter density
    fprintf(fid, '%d 1 %f %f %f 1 1 1 \n', ID(i), X(i), Y(i), 0);
end

fprintf(fid, '\nBonds\n\n');
for i = 1:length(idx)
    if ~idx(i)
        counter = counter + 1;
        % Format: bondID  bondType  atom1  atom2
        fprintf(fid, '%d %d %d %d\n', counter, bondType(i), iatom(i), jatom(i));
    end
end

fclose(fid);
fprintf('Wrote %s with %d atoms and %d bonds.\n', lammps_visual_file, Atom_count, Bond_count);

end