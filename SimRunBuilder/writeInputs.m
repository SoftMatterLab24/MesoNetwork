function writeInputs(Data,logs,options)
% writeInputs - Write LAMMPS input files based on template and log data
%
% Syntax: writeInputs(Data,logs,options)
%
% Inputs:
%    Data - Structure containing .dat, .log, and .table file information
%    logs - Structure array containing parsed log file data
%    options - Structure containing simulation options
% Outputs:
%    None (input files are written to disk)
% Example: writeInputs(Data, logs, options)
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

    % A. Grab the correct Input Template
    if strcmpi(options.dist_type,'bimodal')
        if strcmpi(options.sim_type,'unnotched')
            inputLoc = fullfile('../InputScripts/Templates/Bimodal/BD_LoadUnnotch_TEMPLATE.in');
            %editlines = [16 56 61 62 72 76 77];
        elseif strcmpi(options.sim_type,'notched')
            inputLoc = fullfile('../InputScripts/Templates/Bimodal/BD_LoadUnnotch_TEMPLATE.in');
        end
    elseif strcmpi(options.dist_type,'monodisperse')

    elseif strcmpi(options.dist_type,'polydisperse')
    else
    end

    % B. Read necessary files
    %%% ---- Input Template ---- %%%
    fid = fopen(inputLoc,'r');
    inputTemplate = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);

    %%% --- Batch Sript --- %%%
    if options.iHPC
        if strcmpi(options.HPC_type,'Dane')
            batchLoc = fullfile('../InputScripts/Templates/Batch/Dane/DaneBatch_TEMPLATE.pbs');
        elseif strcmpi(options.HPC_type,'Milan')
            batchLoc = fullfile('../InputScripts/Templates/Batch/Milan/MilanBatch_TEMPLATE.pbs');
        else
            error('Unknown HPC_type %s', options.HPC_type)
        end

        fidB = fopen(batchLoc,'r');
        batchTemplate = textscan(fidB, '%s', 'Delimiter', '\n');
        fclose(fidB);
    end

    % C. Loop over logs and write input files
    for iLog = 1:Data.Nlogs

        % -1. Compute properties for clamps
        Domain = logs(iLog).domain_size;
        ylo = Domain(2,1); yhi = Domain(2,2);

        upperClampBound = yhi - options.clampFrac*yhi;
        lowerClampBound = ylo - options.clampFrac*ylo;
        
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

        % 1. Get both bond.table and LD.table - check that sample number and replicate number both match
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

        % 2. Create directory to place input file
        simDir = fullfile(options.write_location, sprintf('%04d', iLog));
        if ~exist(simDir, 'dir')
            mkdir(simDir);
        end

        % 3. Open new input file for writing
        inputSplit = split(inputLoc,'\');
        inputType = inputSplit{end};
        inputParts = split(inputType,'_');

        inputFileName = sprintf('%s_%s_SMP%04d_N%04d.in',inputParts{1},inputParts{2},logs(iLog).sample_number,logs(iLog).replicate_number);
        inputFilePath = fullfile(simDir, inputFileName);
        inFileName{iLog} = inputFilePath;
        fidOut = fopen(inputFilePath, 'w');

        % 4. Modify template lines based on log data
        inputLines = inputTemplate{1};
        if strcmpi(options.dist_type,'bimodal')
            if strcmpi(options.sim_type,'unnotched')
                %%% ----- Unnotched Bimodal ----- %%%
                inputLines{16} = sprintf('read_data           %s', datName);
                inputLines{56} = sprintf('variable            sig equal $(%2.2f*v_b)    #LJ particle diameter',logs(iLog).rc);
                inputLines{61} = sprintf('region              top_clamped block INF INF %4.2f INF INF INF',upperClampBound);
                inputLines{62} = sprintf('region              bot_clamped block INF INF INF %4.2f INF INF',lowerClampBound);
                inputLines{72} = sprintf('pair_coeff          * * %s',LDtab);
                inputLines{76} = sprintf('bond_coeff          1 $(v_kBT) $(100*v_kBT) 0 %s KEY stretch $(v_lamc)',bondtab);
                inputLines{77} = sprintf('bond_coeff          1 $(v_kBT) $(100*v_kBT) 0 %s KEY stretch $(v_lamc)',bondtab);
            elseif strcmpi(options.sim_type,'notched')
                %%% ----- Notched Bimodal ----- %%%
            end
        elseif strcmpi(options.dist_type,'monodisperse')
            if strcmpi(options.sim_type,'unnotched')
                %%% ----- Unnotched Mono ----- %%%
            elseif strcmpi(options.sim_type,'notched')
                %%% ----- Notched Mono ----- %%%
            end
        elseif strcmpi(options.dist_type,'polydisperse')
            if strcmpi(options.sim_type,'unnotched')
                %%% ----- Unnotched Mono ----- %%%
            elseif strcmpi(options.sim_type,'notched')
                %%% ----- Notched Mono ----- %%%
            end
        else
        end

        % 5. Write modified lines to new input file
        for j = 1:length(inputLines)
            fprintf(fidOut, '%s\n', inputLines{j});
        end

        fclose(fidOut); % Close the output file

        % 6. If HPC, write batch script
        if options.iHPC
            batchFileName = fullfile(simDir, sprintf('Submit_%04d.sh', iLog));
            fidBatchOut = fopen(batchFileName, 'w');
            batchLines = batchTemplate{1};
            
            % Modify batch lines depending on type
            if strcmpi(options.HPC_type,'Dane')
                splitBatchLines = split(batchLines{7},' ');
                RunLocPrefix = splitBatchLines{end};            %Split to get path prefix
                RunLoc = fullfile(RunLocPrefix,inputFilePath);  %Create full path to input file
                RunLoc = strrep(RunLoc,'\','/');                %Convert \ to / <- For UNIX systems
                splitBatchLines{end} = RunLoc;                  %Replace last entry with full path

                batchLines{7} = sprintf('%s',string(join(splitBatchLines,' ')));
            elseif strcmpi(options.HPC_type,'Milan')
                splitBatchLines = split(batchLines{17},' ');
                RunLocPrefix = splitBatchLines{end};
                RunLoc = fullfile(RunLocPrefix,inputFilePath);
                RunLoc = strrep(RunLoc,'\','/');                %<- For UNIX systems
                splitBatchLines{end} = RunLoc;

                batchLines{17} = sprintf('%s',string(join(splitBatchLines,' ')));
            else 
                error('Unknown HPC_type %s', options.HPC_type)
            end

            for j = 1:length(batchLines)
                fprintf(fidBatchOut, '%s\n', batchLines{j});
            end
            fclose(fidBatchOut);
        end

        % 7. Copy necessary .table, .dat, and .log files to simDir
        try
            copyfile(fullfile(options.data_location, bondtab), fullfile(simDir, bondtab));
            copyfile(fullfile(options.data_location, LDtab), fullfile(simDir, LDtab));
            copyfile(fullfile(options.data_location, datName), fullfile(simDir, datName));
            copyfile(fullfile(options.data_location, Data.logs{iLog}), fullfile(simDir, Data.logs{iLog}));
            fprintf('Copied files for %s to %s\n', inputFilePath, simDir);  
        catch
            error('Failed to copy files')
        end
        
    end

    % D. Write bash script to queue all jobs (if HPC)
    if options.iHPC
        bashFileName = fullfile(options.write_location, 'BigBrother.sh');
        fidBashOut = fopen(bashFileName, 'w');
        fprintf(fidBashOut, '#!/bin/bash\n\n');
        for iLog = 1:Data.Nlogs
            simDir = fullfile(options.write_location, sprintf('%04d', iLog));
            batchFileName = fullfile(simDir, sprintf('Submit_%04d.sh', iLog));
            fprintf(fidBashOut, 'cd %s\n', simDir);
            fprintf(fidBashOut, 'sbatch %s\n', batchFileName);
            fprintf(fidBashOut, 'cd ..\n\n');
        end
        fclose(fidBashOut);
        fprintf('Wrote batch submission script: %s\n', bashFileName);
    end
end