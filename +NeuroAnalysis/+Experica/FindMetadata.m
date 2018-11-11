function [i] = FindMetadata(metatable,test,range)
%FINDMETADATA Search for matching tests in range, return indices
%   Detailed explanation goes here

%% Try Find Optical Imaging First
i = metatable.iquery({'sourceformat'}, {'OI'},range);
if ~isempty(i)
    return;
end
%% No Optical Imaging Found, Find Experica
searchtemplate = NeuroAnalysis.Base.getstructfields(test,...
    {'ID', 'RecordSite','sourceformat','filename'});
keys = fieldnames(searchtemplate);
values = struct2cell(searchtemplate);
i = metatable.iquery(keys, values,range);

end