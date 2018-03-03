function [meta] = PrepareMetadata( dataset,callbackresult )
%PREPAREMETADATA Prepare Optical Imaging metadata
%   Detailed explanation goes here

meta=[];
if (isempty(dataset) || ~isfield(dataset,'image'))
    return;
end

fields = {'Subject_ID','RecordSession'};
meta = NeuroAnalysis.Base.getstructfields(dataset,fields);

meta.files={dataset.filepath};
meta.sourceformat = dataset.sourceformat;

end

