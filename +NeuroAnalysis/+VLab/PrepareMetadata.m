function [meta] = PrepareMetadata( dataset,callbackresult )
%PREPAREMETADATA Prepare VLab metadata
%   Detailed explanation goes here

meta=[];
if (isempty(dataset) || ~isfield(dataset,'ex'))
    return;
end

fields = {'ID', 'Subject_ID','RecordSite','RecordSession'};
meta = NeuroAnalysis.Base.getstructfields(dataset.ex,fields);

end

