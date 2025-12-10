function writeInputs(Data,logs,options)
    
    % 0. Create directory to place 

    % A. Grab the correct Input Template
    if strcmpi(options.dist_type,'bimodal')
        if strcmpi(options.sim_type,'unnotched')
            inputLoc = fullfile('../InputScripts/Templates/Bimodal/BD_LoadUnnotch_TEMPLATE.in');
            editlines = [16 56 61 62 72 76 77];
        elseif strcmpi(options.sim_type,'notched')
        end
    elseif strcmpi(options.dist_type,'monodisperse')

    elseif strcmpi(options.dist_type,'polydisperse')
    else


    end

    % B. Read the template file
    fid = fopen(inputLoc,'r');
    inputTemplate = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);

    % C. Loop over logs and write input files
    for iLog = 1:Data.Nlogs
        
        % 0. Get the name of .dat file from Data structure - check that sample number and replicate number both match
        datName = [];
        for idat = 1:Data.Ndats
            % parse sample and replicate from dat file name
            datNameParts = split(Data.dats{idat},'_');
            sampleNum_str = datNameParts{end-1};
            replicateNum_str = datNameParts{end};
            sampleNum = str2double(sampleNum_str(4:end));
            replicateNum = str2double(replicateNum_str(2:end-4));
            if sampleNum == logs(iLog).sample_number && replicateNum == logs(iLog).replicate_number
                datName = Data.dats{idat};
                break;
            end
        end

        % Get both bond.table and LD.table - check that sample number and replicate number both match
        bondtab = []; LDtab = [];
        for idat = 1:Data.Ntabs
            % parse sample and replicate from table file name
            tabNameParts = split(Data.tables{idat},'_');
            sampleNum_str = tabNameParts{end-1};
            replicateNum_str = tabNameParts{end};
            sampleNum = str2double(sampleNum_str(4:end));
            replicateNum = str2double(replicateNum_str(2:5));
            if sampleNum == logs(iLog).sample_number && replicateNum == logs(iLog).replicate_number
                if contains(Data.tables{idat},'bond')
                    bondtab = Data.tables{idat};
                elseif contains(Data.tables{idat},'LD')
                    LDtab = Data.tables{idat};
                end
            end
        end

        % 1. Create directory to place input file
        simDir = fullfile(options.write_location, sprintf('%04d', iLog));
        if ~exist(simDir, 'dir')
            mkdir(simDir);
        end

        % 2. Open new input file for writing
        inputFileName = fullfile(simDir, sprintf('Input_Sim_%d.in', iLog));
        inFileName{iLog} = inputFileName;
        fidOut = fopen(inputFileName, 'w');

        % 3. Modify template lines based on log data
        inputLines = inputTemplate{1};
        if strcmpi(options.dist_type,'bimodal')
            %%% ----- Unnotched Bimodal ----- %%%
            if strcmpi(options.sim_type,'unnotched')
                inputLines{16} = sprintf('read_data           %s', datName);
                inputLines{56} = sprintf('variable            sig equal $(%2.2f*v_b)    #LJ particle diameter',logs(iLog).rc);
                % finish editing other lines
                %editlines = [16 56 61 62 72 76 77];
            %%% ----- Notched Bimodal ----- %%%
            elseif strcmpi(options.sim_type,'notched')
            end
        elseif strcmpi(options.dist_type,'monodisperse')
        elseif strcmpi(options.dist_type,'polydisperse')
        else
        end

        % 4. Write modified lines to new input file
        for j = 1:length(inputLines)
            fprintf(fidOut, '%s\n', inputLines{j});
        end

        % 5. Close the output file
        fclose(fidOut);

        % Copy necessary .table, .dat, and .log files to simDir
        try
            copyfile(fullfile(options.data_location, bondtab), fullfile(simDir, bondtab));
            copyfile(fullfile(options.data_location, LDtab), fullfile(simDir, LDtab));
            copyfile(fullfile(options.data_location, datName), fullfile(simDir, datName));
            copyfile(fullfile(options.data_location, Data.logs{iLog}), fullfile(simDir, Data.logs{iLog}));
            fprintf('Copied files for %s to %s\n', inputFileName, simDir);  
        catch
            error('Failed to copy files')
        end
        
    end

    

end