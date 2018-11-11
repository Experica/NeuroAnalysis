function [ result ] = CollectMetadata( datafile,isparallel,varargin )
%COLLECTMETADATA Collect metadata from exported dataset
%   Detailed explanation goes here

%% Batch collect
if isa(datafile,'cell')
    funlist=repelem({'NeuroAnalysis.IO.CollectMetadata'},length(datafile));
    vararginlist = arrayfun(@(i)[i,{isparallel},varargin],datafile,'UniformOutput',false);
    result = NeuroAnalysis.Base.ApplyFunctions(funlist,vararginlist,isparallel);
    return;
end
%% Open dataset
result.status = false;
result.source = datafile;
dataset=[];
result.meta = [];
if ~isa(datafile,'char') || exist(datafile, 'file') ~= 2
    return;
end
disp(['Reading Dataset:    ',datafile,'    ...']);
dataset = NeuroAnalysis.IO.loaddataset(datafile,{'any data fields not needed for metadata collection.'});
disp(['Reading Dataset:    ',datafile,'    Done.']);
if isempty(dataset)
    return;
end
%% Extract metadata
meta = NeuroAnalysis.Base.EvalFun(...
    ['NeuroAnalysis.',dataset.sourceformat,'.PrepareMetadata'],...
    {dataset,varargin});
if (~isempty(meta))
    result.meta = meta;
    result.status = true;
end

end