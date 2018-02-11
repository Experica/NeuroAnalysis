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
result.source = '';
if (~isa(datafile,'char'))
    return;
end
if ~strcmp(sourceformat,'Unknown')
    % Prepare dataset from file
    dataset = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.',sourceformat,'.Prepare'],...
        [{datafile, exportdir},varargin]);
    if isempty(dataset) || (isfield(dataset,'status') && ~dataset.status)
        return;
    end
    % Get the export path
    if ~isfield(dataset, 'filepath') || isempty(dataset.filepath)
        dataset.filepath = fullfile(exportdir,[NeuroAnalysis.Base.filenamenodirext(datafile),'.mat']);
    end
    exportpath = dataset.filepath;
else
    % Unknown source, just prepare the path for copied file
    [~,filename,ext]=fileparts(datafile);
    exportpath=fullfile(exportdir,[filename,ext]);
end
result.source = datafile;
%% Save dataset
disp(['Exporting Dataset:    ',exportpath,'    ...']);
if ~strcmp(sourceformat,'Unknown')
    save(exportpath,'dataset','-v7.3');
else
    copyfile(datafile,exportpath);
end
disp('Exporting Dataset:    Done.');
%% Callback
callbackfun = callback{1};
callbackarg = callback{2};
callbackresult = struct([]);
if ~isempty(callbackfun) && ~strcmp(sourceformat,'Unknown')
    disp(['Applying Callback:    ',callbackfun,'    ...']);
    callbackresult = NeuroAnalysis.Base.EvalFun(callbackfun,[{dataset},callbackarg]);
    disp('Applying Callback:    Done.');
end
%% Prepare Metadata
result.meta = struct([]);
if ~strcmp(sourceformat,'Unknown')
    metaresult = NeuroAnalysis.Base.EvalFun( ...
        ['NeuroAnalysis.',sourceformat,'.PrepareMetadata'], ...
        {dataset,callbackresult});
    if metaresult.status
        result.meta = metaresult.meta;
    end
end
result.status = true;
