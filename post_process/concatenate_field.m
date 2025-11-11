function out = concatenate_field(S, fieldname)
% Collects S(:).fieldname vertically
out = [];
for k=1:numel(S), out = [out; S(k).(fieldname)]; end %#ok<AGROW>
end