function order = ComputeOrder(obj, Atoms, Bonds)
% -------------------------------------------------------------------------
% ComputeOrder
% - Compute structural order parameters for the generated network.
%
% Atoms column layout: [ID | molID | X | Y | Z | deg | nbrs...]
% X at col 3, Y at col 4.
% -------------------------------------------------------------------------

    %#ok<INUSD>
    order = struct();

    [phi6k, phi6_hexatic, phi6_hexagonal] = ComputeHexOrder(Atoms, Bonds);

    order.hex.phi6k           = phi6k;
    order.hex.phi6_hexatic    = phi6_hexatic;
    order.hex.phi6_hexagonal  = phi6_hexagonal;

    obj.log.print('   Computed structural order parameters:\n');
    obj.log.print('   Hexatic order phi6 = %.4f\n', phi6_hexatic);
    obj.log.print('   Hexagonal order phi6 = %.4f\n', phi6_hexagonal);

end


function [phi6k, phi6_hexatic, phi6_hexagonal] = ComputeHexOrder(Atoms, Bonds)
% Compute local and global hexagonal/hexatic order parameters from
% network connectivity. Uses X/Y at cols 3/4 under the new Atoms layout.

    N = size(Atoms,1);

    phi6k    = zeros(N,1);
    phi6norm = zeros(N,1);

    if N == 0
        phi6_hexatic = 0;
        phi6_hexagonal = 0;
        return;
    end

    neighborList = cell(N,1);

    for iB = 1:size(Bonds,1)
        a1 = Bonds(iB,2);
        a2 = Bonds(iB,3);

        if a1 >= 1 && a1 <= N && a2 >= 1 && a2 <= N
            neighborList{a1} = [neighborList{a1}, a2]; %#ok<AGROW>
            neighborList{a2} = [neighborList{a2}, a1]; %#ok<AGROW>
        end
    end

    for k = 1:N

        neighbors = neighborList{k};
        numNeighbors = numel(neighbors);

        if numNeighbors == 0
            phi6k(k) = 0;
            phi6norm(k) = 0;
            continue;
        end

        angles = zeros(numNeighbors,1);

        for j = 1:numNeighbors
            vec = Atoms(neighbors(j), 3:4) - Atoms(k, 3:4);   % X,Y at cols 3,4
            angles(j) = atan2(vec(2), vec(1));
        end

        phi6k(k) = sum(exp(1i * 6 * angles)) / numNeighbors;
        phi6norm(k) = abs(phi6k(k));
    end

    phi6_hexatic   = abs(sum(phi6k) / N);
    phi6_hexagonal = sum(phi6norm) / N;

end
