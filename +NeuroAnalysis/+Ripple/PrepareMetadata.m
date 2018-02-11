function [ result ] = PrepareMetadata( dataset,callbackresult )
%PREPAREMETADATA Prepare metadata for exported dataset and its callbackresult
%   Detailed explanation goes here

result=[];
%% Prepare experimental metadata
exmeta = [];
if (~isempty(dataset) && isfield(dataset,'ex'))
    exmeta = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.',dataset.ex.sourceformat,'.PrepareMetadata'], ...
        {dataset,callbackresult});
end
%% Prepare Ripple metadata
ripplemeta = [];
if (~isempty(dataset))
    ripplemeta.files={dataset.filepath};
    [~, ripplemeta.filename, ~] = fileparts(dataset.filepath);
    ripplemeta.sourceformat = dataset.sourceformat;
end
%% Get Callback metadata
callbackmeta = [];
if(~isempty(callbackresult) && isfield(callbackresult,'result'))
    callbackmeta = callbackresult.result;
end
%% Combine metadata
meta = NeuroAnalysis.Base.mergestruct(exmeta,ripplemeta);
if (~isempty(meta) && ~isempty(fieldnames(meta)))
    if (~isempty(callbackmeta) && ~isempty(fieldnames(callbackmeta)))
        meta.result = callbackmeta;
    end
    result=meta;
end

end

