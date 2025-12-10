function [Data] = ParseDirectory(data_location)
    
listing = dir(data_location);

% Find all .log files in the directory

indDat = 0; indLog = 0; indTab = 0;
for ii = 1:length(listing)
    name = listing(ii).name;
    
    %find .log files
    if contains(name,'.log')
        indLog = indLog + 1;
        logs{indLog} = name;
        continue
    end
    
    %find .table files
    if contains(name,'.table')
        indTab = indTab + 1;
        tables{indTab} = name;
        continue
    end

    %find .dat files 
    if contains(name,'.dat')
        indDat = indDat +1;
        dats{indDat} = name;
        continue
    end
end

% Write parsed info to data structure
Data.Nlogs = indLog;
Data.Ndats = indDat;
Data.Ntabs = indTab;
Data.logs = logs;
Data.tables = tables;
Data.dats = dats;

end