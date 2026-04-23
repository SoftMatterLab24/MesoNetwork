function [Atoms, LatticeData] = AddAtomsHex(obj)
% -------------------------------------------------------------------------
% AddAtomsHex
% - Generate a 2D hexagonal (triangular) lattice
% - No geometric disorder applied here
%
% OUTPUT:
%   Atoms       : [ID | molID | X | Y | Z | deg | nbr_1 ... nbr_{Max_peratom_bond}]
%                 All atoms receive molID = 1 (single-molecule default).
%   LatticeData : struct with idx_map, Nx, Ny
% -------------------------------------------------------------------------

    xlo = obj.domain.xlo;
    xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;
    yhi = obj.domain.yhi;

    Lx = xhi - xlo;
    Ly = yhi - ylo;

    a = obj.domain.min_node_sep;

    dy = a * sqrt(3) / 2;

    Ny_est = floor(Ly / dy) + 2;
    Nx_est = floor(Lx / a) + 3;

    idx_map = zeros(Ny_est, Nx_est);

    maxNodes = Ny_est * Nx_est;
    x_all = zeros(maxNodes,1);
    y_all = zeros(maxNodes,1);

    nat = 0;

    for iy = 1:Ny_est

        y = ylo + (iy-1) * dy;

        if (y < ylo) || (y > yhi)
            continue;
        end

        x_offset = 0;
        if mod(iy,2) == 1
            x_offset = a/2;
        end

        for ix = 1:Nx_est

            x = xlo + (ix-1)*a + x_offset;

            if (x < xlo) || (x > xhi)
                continue;
            end

            nat = nat + 1;
            x_all(nat) = x;
            y_all(nat) = y;
            idx_map(iy, ix) = nat;

        end
    end

    x_all = x_all(1:nat);
    y_all = y_all(1:nat);

    Max_peratom_bond = obj.peratom.Max_peratom_bond;
    ncols = 6 + Max_peratom_bond;

    Atoms = zeros(nat, ncols);
    Atoms(:,1) = (1:nat).';   % ID
    Atoms(:,2) = 1;            % molID (single-molecule default)
    Atoms(:,3) = x_all;        % X
    Atoms(:,4) = y_all;        % Y
    Atoms(:,5) = 0.0;          % Z
    % col 6 = degree (0), cols 7..6+Max_peratom_bond = neighbors (0)

    LatticeData.idx_map = idx_map;
    LatticeData.Nx      = size(idx_map,2);
    LatticeData.Ny      = size(idx_map,1);

end
