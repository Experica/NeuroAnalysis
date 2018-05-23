function [ result ] = Export(datafile, exportdir,sourceformat,isparallel,callback,varargin )
%EXPORT Export prepared dataset in Matlab MAT format file
%   Detailed explanation goes here

%% Batch export
if isa(datafile,'cell')
    funlist=repelem({'NeuroAnalysis.IO.Export'},length(datafile));
    vararginlist = arrayfun(@(i)[i,{exportdir,sourceformat,isparallel,callback},varargin],datafile,'UniformOutput',false);
    result = NeuroAnalysis.Base.ApplyFunctions(funlist,vararginlist,isparallel);
    return;
end
%% Prepare
result.status = false;
result.source = datafile;
dataset=[];
if (~isa(datafile,'char'))
    return;
end
if ~strcmp(sourceformat,'Unknown')
    % Prepare dataset from file
    dataset = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.',sourceformat,'.Prepare'],...
        [{datafile},varargin]);
    if isempty(dataset) || (isfield(dataset,'status') && ~dataset.status)
        return;
    end
    % Get the export path
    if ~isfield(dataset, 'filepath') || isempty(dataset.filepath)
        dataset.filepath = fullfile(exportdir,[NeuroAnalysis.Base.filenamenodirext(datafile),'.mat']);
    end
    exportpath = dataset.filepath;
else
    % Unknown source, just prepare the path for the file to be copied
    [~,filename,ext]=fileparts(datafile);
    exportpath=fullfile(exportdir,[filename,ext]);
end
%% Save dataset
disp(['Exporting Dataset:    ',exportpath,'    ...']);
if ~strcmp(sourceformat,'Unknown')
    save(exportpath,'-struct','dataset','-v7.3');
else
    copyfile(datafile,exportpath);
end
disp('Exporting Dataset:    Done.');
%% Callback
callbackfun = callback{1};
callbackarg = callback{2};
callbackresult = [];
if ~isempty(callbackfun) && ~strcmp(sourceformat,'Unknown')
    disp(['Applying Callback:    ',callbackfun,'    ...']);
    callbackresult = NeuroAnalysis.Base.EvalFun(callbackfun,[{dataset},callbackarg]);
    disp('Applying Callback:    Done.');
end
%% Prepare Metadata
result.meta =[];
if ~strcmp(sourceformat,'Unknown')
    meta = NeuroAnalysis.Base.EvalFun(...
        ['NeuroAnalysis.',sourceformat,'.PrepareMetadata'], ...
        {dataset,callbackresult});
    if (~isempty(meta))
        result.meta = meta;
    end
end
result.status = true;
