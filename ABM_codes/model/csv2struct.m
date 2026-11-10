function T = csv2struct(fn)
%CSV2STRUCT  Read a CSV file written by struct2csv into a struct array.
%   Numeric columns become doubles; other columns stay as character arrays.
fid = fopen(fn, 'r');
if fid < 0, error('Cannot open %s', fn); end
hdr = strtrim(fgetl(fid)); f = strsplit(hdr, ',');
L = {};
while true
    s = fgetl(fid); if ~ischar(s), break; end
    if ~isempty(strtrim(s)), L{end+1} = s; end %#ok<AGROW>
end
fclose(fid);
n = numel(L); C = cell(n, numel(f));
for i = 1:n
    p = strsplit(L{i}, ',');
    for j = 1:numel(f)
        v = str2double(p{j});
        if isnan(v) && ~strcmpi(strtrim(p{j}), 'NaN'), C{i,j} = p{j}; else, C{i,j} = v; end
    end
end
T = cell2struct(C, f, 2);
end
