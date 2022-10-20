function [meta] = readmeta(filepath)
% Parse ini file into cell entries C{1}{i} = C{2}{i}
fid = fopen(filepath, 'r');
C = textscan(fid, '%[^=] = %[^\r\n]');
fclose(fid);

meta = struct();
% Convert each cell entry into a struct entry
for tagi = 1:length(C{1})
    tag = C{1}{tagi};
    if tag(1) == '~'
        % remake tag excluding first character
        tag = sprintf('%s', tag(2:end));
    end
    v = C{2}{tagi};
    [t,r]=str2num(v);
    if r
        v=t;
    end
    meta.(tag) = v;
end
t = meta.fileSizeBytes/meta.nSavedChans;
if mod(t,2) ~= 0
    warning('Binary file does not have same number of samples for all channels. Use the maximum same length for all channels.');
    meta.nFileSamp = int64(floor(t/2));
else
    meta.nFileSamp = int64(t/2);
end
meta.nSavedChans = int64(meta.nSavedChans);
meta.fileDate = datenum(meta.fileCreateTime,'yyyy-mm-ddTHH:MM:SS');

end

