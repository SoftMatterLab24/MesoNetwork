classdef network < handle

properties

    % Define properites and set default parameters
    %options = struct();

    %%% Flags
    flags = struct(...
        'isave', true, ...
        'imanualseed', false ...
    );

    %%% Domain
    domain = struct( ...
        );

    %%% Architecture
 

    %%% Per/atom

    %%% Per/bond

    %%% Defect

    %%% Potential


end

methods

    function [] = generateNetwork(obj)

        %%% Loop over replicates
        for ii=1 %Nreps

            % -----------------------------------%
            %1. Prepare replicate-specific information

            %2. Constuct domain

            %3. Add atoms

            %4. Assign per/atom

            %5. Add bonds

            %6. Assign per/bond

            %7. Add defects

            %8. Clean-up network


            % -----------------------------------%
            %9. Construct local density potential

            %10. Scale domain if needed
            
            % -----------------------------------%
            %11. Show visualization and statistics

            %12. Computes

            %13. Write data files

            % -----------------------------------%

        end

    end


end
end