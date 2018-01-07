function [ result ] = Export(dataset, exportdir,sourceformat,isparallel,callback,varargin )
%EXPORT Export prepared dataset in Matlab MAT format file
%   Detailed explanation goes here

%% Batch export
if isa(dataset,'cell')
    funlist=repelem({'NeuroAnalysis.IO.Export'},length(dataset));
    vararginlist = arrayfun(@(i)[i,{exportdir,sourceformat,isparallel,callback},varargin],dataset,'UniformOutput',false);
    result = NeuroAnalysis.Base.ApplyFunctions(funlist,vararginlist,isparallel);
    return;
end
%% Get dataset and export path
result.status = false;
result.source = '';
if isa(dataset,'char')
    result.source = dataset;
    if ~strcmp(sourceformat,'Unknown')
        dataset = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.',sourceformat,'.Prepare'],[{dataset},varargin]);
        if ~isempty(dataset) && isfield(dataset,'status') && ~dataset.status
            return;
        end
        filename = NeuroAnalysis.Base.filenamenodirext(dataset.source);
        exportpath = fullfile(exportdir,[filename,'.mat']);
    else
        [~,filename,ext]=fileparts(dataset);
        exportpath=fullfile(exportdir,[filename,ext]);
    end
end
%% Save dataset
disp(['Exporting Dataset:    ',exportpath,'    ...']);
if ~strcmp(sourceformat,'Unknown')
    dataset.filepath = exportpath;
    save(exportpath,'dataset','-v7.3');
else
    copyfile(dataset,exportpath);
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
