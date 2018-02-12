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
result.source = '';
dataset=[];
result.meta = [];
if (~isa(datafile,'char'))
    return;
end
disp(['Reading Dataset:    ',datafile,'    ...']);
loadresult = load(datafile);
disp('Reading Dataset:    Done.');
if (~isempty(loadresult) && isfield(loadresult,'dataset'))
    dataset = loadresult.dataset;
else
    return;
end
result.source = datafile;
%% Extract metadata
meta = NeuroAnalysis.Base.EvalFun(...
    ['NeuroAnalysis.',dataset.sourceformat,'.PrepareMetadata'],...
    {dataset,varargin});
if (~isempty(meta))
    result.meta = meta;
    result.status = true;
end

end
