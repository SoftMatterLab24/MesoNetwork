classdef network < handle

properties

    % Define properites and set default parameters

    %%% Flags
    flags = struct(...
        'isave', true, ...
        'imanualseed', false ...
    );

    %%% Domain
    domain = struct(...
        );

    %%% Architecture subclass
    arch architecture 

    %%% Per/atom
    peratom = struct(...
    );

    %%% Per/bond subclass
    perbond bondstyle

    %%% Defect
    defect = struct(...
    );

    %%% Potential
    pot = struct(...
    );

    %%% log
    log = struct(...
    );

end

methods

    % Explicitly initialize subclass objects in the constructor
    function obj = network()
        obj.arch = architecture();
        obj.perbond = bondstyle();
    end

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