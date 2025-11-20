function [lo, hi] = readBoundsLine(fid)
    lo = NaN; hi = NaN;
    while true
        line = fgetl(fid);
        if ~ischar(line)
            return; % EOF, caller handles it
        end
        vals = sscanf(line, '%f');
        if numel(vals) >= 2
            lo = vals(1);
            hi = vals(2);
            return;
        end
        % if line wasn’t numeric (blank/comment), just keep reading
    end
end
