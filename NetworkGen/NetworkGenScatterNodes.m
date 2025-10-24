function [Atoms] = NetworkGenScatterNodes(Domain)
% Atoms layout:
% [ ID | X | Y | Z | num_bond | nbr1 | nbr2 | nbr3 | nbr4 | spare ]

xlo = Domain.xlo; xhi = Domain.xhi;
ylo = Domain.ylo; yhi = Domain.yhi;
zlo = Domain.zlo; zhi = Domain.zhi;

Max_atom                        = 73335; %Domain.Max_atom;
node_scatter_max_tries          = Domain.node_scatter_max_tries;
min_node_sep                    = Domain.min_node_sep;
Rmax                            = min_node_sep;

%%%
Atoms  = zeros(Max_atom, 10);

tic

icount = 0; N_atom = 0; ng = 1;
while (N_atom < Max_atom)
    
    icount = icount + 1;
    if icount > node_scatter_max_tries
        fprintf("   Placed %d atoms out of %d requested in %d\n", N_atom, Max_atom,toc);
        break;
    end

    %%%% Always place first atom
    if N_atom == 0
        xg = (xhi-xlo).*rand(1,1) + xlo;
        yg = (yhi-ylo).*rand(1,1) + ylo;
        zg = 0;

        N_atom = N_atom + 1
        ID = N_atom;
        Xg = xg;
        Yg = yg;
        Zg = zg;
    end

    %%% Randomly place 200 atoms
    xg = (xhi-xlo).*rand(200,1) + xlo;
    yg = (yhi-ylo).*rand(200,1) + ylo;
    zg = zeros(200,1);

    %%% Remove atoms that are too close to atoms in remaining seed
    Idx = rangesearch([xg yg zg],[xg yg zg],Rmax);
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
    Idx = rangesearch([Xg(1:ng) Yg(1:ng) Zg(1:ng)],[xg yg zg],Rmax);
    idelete = ~(cellfun('isempty',Idx));
    if ~isempty(idelete)
        xg(idelete) = [];
        yg(idelete) = [];
        zg(idelete) = [];
    end

    if isempty(xg)
        continue
    end

    N_atom = N_atom + length(xg)

    %%% Add remaining particles to discretization %%%
    n_add = length(xg);
    id = [(ng+1):(ng+n_add)]';
    ID = [ID; id];
    Xg = [Xg; xg];
    Yg = [Yg; yg];
    Zg = [Zg; zg];
    ng = ng+n_add;

    icount = 0;
end

if N_atom == 0
    error('Error: No atoms were placed in the domain. Consider increasing domain size or decreasing minimum distance between atoms.');
end

    % Write to atom vector
    Atoms(1:N_atom,1) = ID;
    Atoms(1:N_atom,2) = Xg;
    Atoms(1:N_atom,3) = Yg;
    Atoms(1:N_atom,4) = Zg;

end