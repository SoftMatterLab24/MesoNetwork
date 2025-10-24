function [coords,VolGrain] = fill_stl_fun(loc,file,Rmean,scale)

%% Open .stl file

tic

%% Need to do a dimensions check sometimes (mm or um)

[~,coords] = stlread(loc+'\'+file);  % stlread sometimes doesn't work...
coords     = scale*coords; %if mm need to convert to um
cp0        = mean(coords);
coords     = coords-cp0;

%% Fill sand grain with discrete particles

%%% Scale hull of grain down one particle radius %%%
scale = 1.0;
means = 0;
while means < Rmean
    scale = scale-0.0001;
    means = mean(sqrt(sum((coords - (coords*scale)).^2,2)));
    if scale < 0
        error('Bad scaling')
    end
end
coords = coords*scale;

%%% Create convex hull and alphaShape %%%
xv  = coords(:,1);
yv  = coords(:,2);
zv  = coords(:,3);
coordsv = unique([xv yv zv],'rows','stable');
xv  = coordsv(:,1);
yv  = coordsv(:,2);
zv  = coordsv(:,3);
shp = alphaShape(xv,yv,zv,Inf);

%%% Calculate dimensions of convex hull
Vmean    = (4/3)*pi*Rmean^3;
VolGrain = volume(shp);
Ntry     = floor(VolGrain/Vmean);
Ntry     = 1e5;

xmin = min(xv);
xmax = max(xv);

ymin = min(yv);
ymax = max(yv);

zmin = min(zv);
zmax = max(zv);

%% Fill particles

Rvo2  = xv.^2 + yv.^2 + zv.^2;
Rmax2 = max(Rvo2);

Xg = zeros(Ntry,1);
Yg = zeros(Ntry,1);
Zg = zeros(Ntry,1);

%%% Place first point randomly inside convex hull %%%
in  = 0;
while ~in
    xg = (xmax-xmin).*rand(1,1) + xmin;
    yg = (ymax-ymin).*rand(1,1) + ymin;
    zg = (zmax-zmin).*rand(1,1) + zmin;
    in = inShape(shp,xg,yg,zg);
end
Xg(1) = xg;
Yg(1) = yg;
Zg(1) = zg;

Dmax = (2.05*Rmean);
ng = 1;
icount = 0;
while ng < Ntry

    icount = icount+1;
    if icount > 50
        fprintf('Filled %d particles in %g seconds\n',ng,toc)
        break
    end

    %%% Randomly insert 200 particles %%%
    xg = (xmax-xmin).*rand(200,1) + xmin;
    yg = (ymax-ymin).*rand(200,1) + ymin;
    zg = (zmax-zmin).*rand(200,1) + zmin;

    %%% Remove particles outside of bounding sphere %%%
    Rg2 = xg.^2 + yg.^2 + zg.^2; 
    in  = Rg2 > Rmax2; 
    xg(in) = [];
    yg(in) = [];
    zg(in) = [];

    %%% Remove particles outside of grain convex hull %%%
    in = inShape(shp,xg,yg,zg);
    xg(~in) = [];
    yg(~in) = [];
    zg(~in) = [];

    if isempty(xg)
        continue
    end

    %%% Remove overlapping particles in remaining seeding %%%
    Idx = rangesearch([xg yg zg],[xg yg zg],Dmax);
    for II = 1:length(Idx)
        if Idx{II}(1) == II
            Idx{II}(1) = [];
        end
    end
    i_overlap = ~(cellfun('length',Idx)) == 0;
    if sum(i_overlap)~=0
        Idx = Idx(i_overlap);
        Idx = [Idx{:}]';
        if ~isempty(Idx)
            xg(Idx) = [];
            yg(Idx) = [];
            zg(Idx) = [];
        end
    end

    if isempty(xg)
        continue
    end

    %%% Remove particles overlapping with current discretization %%%
    Idx = rangesearch([Xg(1:ng) Yg(1:ng) Zg(1:ng)],[xg yg zg],Dmax);
    idelete = ~(cellfun('isempty',Idx));
    if ~isempty(idelete)
        xg(idelete) = [];
        yg(idelete) = [];
        zg(idelete) = [];
    end

    if isempty(xg)
        continue
    end

    %%% Add remaining particles to discretization %%%
    n_add = length(xg);
    Xg((ng+1):(ng+n_add)) = xg;
    Yg((ng+1):(ng+n_add)) = yg;
    Zg((ng+1):(ng+n_add)) = zg;
    ng = ng+n_add;

    icount = 0;
end

coords = [Xg(1:ng) Yg(1:ng) Zg(1:ng)];
coords = coords-mean(coords);
coords = coords+cp0;

end