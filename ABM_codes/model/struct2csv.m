function struct2csv(fn, T)
%STRUCT2CSV  Write a struct array of scalar/char fields to a CSV file.
f = fieldnames(T);
fid = fopen(fn, 'w');
fprintf(fid, '%s', f{1}); fprintf(fid, ',%s', f{2:end}); fprintf(fid, '\n');
for i = 1:numel(T)
    for j = 1:numel(f)
        v = T(i).(f{j});
        if j > 1, fprintf(fid, ','); end
        if ischar(v), fprintf(fid, '%s', v);
        elseif isempty(v) || ~isfinite(v), fprintf(fid, 'NaN');
        else, fprintf(fid, '%.6g', v); end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
