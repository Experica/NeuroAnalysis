function [ result ] = Export(dataset, exportpath,sourceformat,isparallel,callback,varargin )
%EXPORT Export prepared dataset in Matlab MAT format file
%   Detailed explanation goes here

%% Batch export
if isa(dataset,'cell')
    funlist=repelem({'NeuroAnalysis.IO.Export'},length(dataset));
    vararginlist = arrayfun(@(i)[i,{exportpath,sourceformat,isparallel,callback},varargin],dataset,'UniformOutput',false);
    result = NeuroAnalysis.Base.ApplyFunctions(funlist,vararginlist,isparallel);
    return;
end
%% Get dataset and export path
if isa(dataset,'char')
    filename = dataset;
    if ~strcmp(sourceformat,'Unknown')
        dataset = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.',sourceformat,'.Prepare'],[{dataset},varargin]);
        if ~isempty(dataset) && isfield(dataset,'status') && ~dataset.status
            result = dataset;
            result.source = NeuroAnalysis.Base.filenamenodirext(filename);
            return;
        end
        filename = NeuroAnalysis.Base.filenamenodirext(dataset.source);
        filepath = fullfile(exportpath,[filename,'.mat']);
    else
        [~,filename,ext]=fileparts(dataset);
        filepath=fullfile(exportpath,[filename,ext]);
    end
end
%% Save dataset
disp(['Exporting Dataset:    ',filepath,'    ...']);
if ~strcmp(sourceformat,'Unknown')
    dataset.filepath = filepath;
    save(filepath,'dataset','-v7.3');
else
    copyfile(dataset,filepath);
    dataset = struct;
    dataset.filepath = filepath;
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
%% Update metadata file
result = NeuroAnalysis.Base.EvalFun('NeuroAnalysis.IO.MetadataExtract', {dataset,exportpath,sourceformat,isparallel,callbackresult});
result.status = true;
result.source = dataset.filepath;
