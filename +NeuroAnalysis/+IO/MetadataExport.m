function [ result ] = MetadataExport( filepath,isparallel,varargin )
%METADATAUPDATE Update metadata for exported files
%   Detailed explanation goes here

%% Batch update
if isa(filepath,'cell')
    funlist=repelem({'NeuroAnalysis.IO.MetadataExport'},length(filepath));
    vararginlist = arrayfun(@(i)[i,{isparallel},varargin],filepath,'UniformOutput',false);
    result = NeuroAnalysis.Base.ApplyFunctions(funlist,vararginlist,isparallel);
    return;
end

%% Open dataset
result.status = false;
result.source = '';
result.meta = struct([]);
assert(isa(filepath,'char'))
disp(['Reading Dataset:    ',filepath,'    ...']);
result.source = filepath;
try
    datafile = load(filepath);
    dataset = datafile.dataset;
catch e
    return % no metadata
end
disp('Reading Dataset:    Done.');
%% Extract metadata
metaresult = NeuroAnalysis.Base.EvalFun(...
    ['NeuroAnalysis.',dataset.sourceformat,'.MetadataPrepare'],...
    {dataset,varargin});
if metaresult.status
    result.meta = metaresult.meta;
    result.status = true;
end
end