function [ meta ] = PrepareMetadata( dataset,callbackresult )
%PREPAREMETADATA Prepare metadata for exported dataset and its callbackresult
%   Detailed explanation goes here

meta=[];
%% Prepare experimental metadata
exmeta = [];
if (~isempty(dataset) && isfield(dataset,'ex'))
    exmeta = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.',dataset.ex.sourceformat,'.PrepareMetadata'], ...
        {dataset,callbackresult});
    if isfield(exmeta,'status') && ~exmeta.status
        exmeta = [];
    end
    if isfield(exmeta,'sourceformat')
        exmeta = rmfield(exmeta,'sourceformat');
    end
    if isfield(exmeta,'filename')
        exmeta = rmfield(exmeta,'filename');
    end
    if isfield(exmeta,'files')
        exmeta = rmfield(exmeta,'files');
    end
end
%% Prepare Ripple metadata
disp('Preparing Ripple Metadata:    ...');
ripplemeta = [];
if (~isempty(dataset))
    ripplemeta.files={dataset.filepath};
    [~, ripplemeta.filename, ~] = fileparts(dataset.filepath);
    ripplemeta.sourceformat = dataset.sourceformat;
end
%% Get Callback metadata
callbackmeta = [];
if(~isempty(callbackresult) && isfield(callbackresult,'result')) && isstruct(callbackresult.result)
    callbackmeta = callbackresult.result;
end
%% Combine metadata
meta = NeuroAnalysis.Base.mergestruct(exmeta,ripplemeta);
if (~isempty(meta) && ~isempty(fieldnames(meta)))
    if (~isempty(callbackmeta) && ~isempty(fieldnames(callbackmeta)))
        meta.result = callbackmeta;
    end
end
disp('Preparing Ripple Metadata:    Done.');
end