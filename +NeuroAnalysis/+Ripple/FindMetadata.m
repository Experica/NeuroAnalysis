function [i] = FindMetadata(metatable,test,range)
%FINDMETADATA Search for matching tests in range, return indices
%   Detailed explanation goes here

searchtemplate = NeuroAnalysis.Base.getstructfields(test,...
    {'ID', 'RecordSite','sourceformat','filename'});
keys = fieldnames(searchtemplate);
values = struct2cell(searchtemplate);
i = metatable.iquery(keys, values,range);

end