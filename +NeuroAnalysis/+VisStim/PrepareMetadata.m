function [meta] = PrepareMetadata( dataset,callbackresult )
%PREPAREMETADATA Prepare VisStim metadata
%   Detailed explanation goes here

meta=[];
if (isempty(dataset) || ~isfield(dataset,'ex'))
    return;
end

fields = {'ID','Subject_ID','RecordSite','RecordSession','File_ID','date'};
meta = NeuroAnalysis.Base.getstructfields(dataset.ex,fields);

end

