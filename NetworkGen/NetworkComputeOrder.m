function [order] = NetworkComputeOrder(Atoms,Bonds)

    % Hexagonal order parameter
    [phi6k,phi6_hexatic,phi6_hexagonal] = ComputeHexOrder(Atoms,Bonds);

    % Other order parameters can be added here

    order.hex.phi6k = phi6k;
    order.hex.phi6_hexatic = phi6_hexatic;
    order.hex.phi6_hexagonal = phi6_hexagonal;
end

function [phi6k,phi6_hexatic,phi6_hexagonal] = ComputeHexOrder(Atoms,Bonds)

    N = size(Atoms,1);
    phi6 = zeros(N,1);
    
    % Build neighbor list
    neighborList = cell(N,1);
    for iB = 1:size(Bonds,1)
        a1 = Bonds(iB,2);
        a2 = Bonds(iB,3);
        neighborList{a1} = [neighborList{a1}, a2];
        neighborList{a2} = [neighborList{a2}, a1];
    end
    
    % Compute hexagonal order parameter for each atom
    for k = 1:N
        neighbors = neighborList{k};
        numNeighbors = length(neighbors);
        %if numNeighbors < 3
        %    phi6(iA) = 0; % Not enough neighbors to define hex order
        %    continue;
        %end
        
        angles = zeros(numNeighbors,1);
        for j = 1:numNeighbors
            vec = Atoms(neighbors(j), 2:4) - Atoms(k, 2:4);
            angles(j) = atan2(vec(2), vec(1));
        end
        
        % Compute the hexagonal order
        phi6k(k) = sum(exp(1i*6*angles))/numNeighbors;
        phi6norm(k) = abs(phi6k(k));
    end

    phi6_hexatic = abs(sum(phi6k)/N);
    phi6_hexagonal = sum(phi6norm)/N;

end