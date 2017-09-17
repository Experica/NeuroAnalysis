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
result.status = false;
result.source = '';
if isa(dataset,'char')
    if ~strcmp(sourceformat,'Unknown')
        dataset = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.',sourceformat,'.Prepare'],[{dataset},varargin]);
        filename = NeuroAnalysis.Base.filenamenodirext(dataset.source);
        exportpath = fullfile(exportpath,[filename,'.mat']);
    else
        [~,filename,ext]=fileparts(dataset);
        exportpath=fullfile(exportpath,[filename,ext]);
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
result.status = true;
result.source = filename;
disp('Exporting Dataset:    Done.');
%% Callback
callbackfun = callback{1};
callbackarg = callback{2};
if ~isempty(callbackfun) && ~strcmp(sourceformat,'Unknown')
    disp(['Applying Callback:    ',callfackfun,'    ...']);
    NeuroAnalysis.Base.EvalFun(callbackfun,[{dataset},callbackarg]);
    disp('Applying Callback:    Done.');
end

