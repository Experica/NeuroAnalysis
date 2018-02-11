function [result] = PrepareMetadata( dataset,callbackresult )
%PREPAREMETADATA Prepare VLab metadata
%   Detailed explanation goes here

result=[];
if (isempty(dataset) || ~isfield(dataset,'ex'))
    return;
end

fields = {'ID', 'Subject_ID','RecordSite','RecordSession'};
result = NeuroAnalysis.Base.getstructfields(dataset.ex,fields);
end

