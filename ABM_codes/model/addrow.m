function T = addrow(T, varargin)
%ADDROW  Append a row (name/value pairs and/or one struct of outputs) to a struct array.
%   T = addrow(T, 'rep', r, 'cfg', 'REF', R) appends labels followed by all fields of R.
row = struct();
k = 1;
while k <= numel(varargin)
    if isstruct(varargin{k})
        f = fieldnames(varargin{k});
        for j = 1:numel(f), row.(f{j}) = varargin{k}.(f{j}); end
        k = k + 1;
    else
        row.(varargin{k}) = varargin{k+1}; k = k + 2;
    end
end
if isempty(T), T = row; else, T(end+1) = orderfields(row, T(1)); end
end
