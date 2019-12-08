function [meta] = PrepareMetadata(dataset,callbackresult)
%PREPAREMETADATA Prepare Scanbox metadata for exported dataset and its callbackresult
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
%% Prepare Scanbox metadata
disp('Preparing Scanbox Metadata:    ...');
Scanboxmeta = [];
if ~isempty(dataset)
    Scanboxmeta.files={dataset.filepath};
    [~, Scanboxmeta.filename, ~] = fileparts(dataset.filepath);
    Scanboxmeta.sourceformat = dataset.sourceformat;
end
%% Get Callback metadata
callbackmeta = [];
if(~isempty(callbackresult) && isfield(callbackresult,'result')) && isstruct(callbackresult.result)
    callbackmeta = callbackresult.result;
end
%% Combine metadata
meta = NeuroAnalysis.Base.mergestruct(exmeta,Scanboxmeta);
if (~isempty(meta) && ~isempty(fieldnames(meta)))
    if (~isempty(callbackmeta) && ~isempty(fieldnames(callbackmeta)))
        meta.result = callbackmeta;
    end
end
disp('Preparing Scanbox Metadata:    Done.');
end

