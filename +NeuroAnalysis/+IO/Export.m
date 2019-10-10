function [ result ] = Export(datafile, exportdir,sourceformat,isparallel,ismergeexport,callback,varargin )
%EXPORT Export prepared dataset in Matlab MAT format file
%   Detailed explanation goes here

global batchexportcallback % does not work in parallel for loop
%% Batch export
if isa(datafile,'cell')
    batchexportcallback=cell(1,2);
    funlist=repelem({'NeuroAnalysis.IO.Export'},length(datafile));
    vararginlist = arrayfun(@(i)[i,{exportdir,sourceformat,isparallel,ismergeexport,callback},varargin],datafile,'UniformOutput',false);
    result = NeuroAnalysis.Base.ApplyFunctions(funlist,vararginlist,isparallel);
    
    callbackfun = batchexportcallback{1};
    if ~isempty(callbackfun)
        disp(['Applying Batch Export Callback:    ',callbackfun,'    ===========================================>']);
        NeuroAnalysis.Base.EvalFun(callbackfun,batchexportcallback(2));
        disp(['Applying Batch Export Callback:    ',callbackfun,'    Done.']);
    end
    clear global batchexportcallback
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
        [{datafile},[varargin,{'exportdir',exportdir}]]);
    if isempty(dataset)
        return;
    else
        if isfield(dataset,'status') && ~dataset.status
            return;
        end
        if isfield(dataset,'earlyfinish') && dataset.earlyfinish
            result.status = true;
            return;
        end
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
    if ismergeexport && exist(exportpath,'file')
        olddataset = matfile(exportpath,'Writable',true);
        fn = fieldnames(dataset);
        for i=1:length(fn)
            olddataset.(fn{i}) = dataset.(fn{i});
        end
    else
        save(exportpath,'-struct','dataset','-v7.3');
    end
else
    copyfile(datafile,exportpath);
end
disp(['Exporting Dataset:    ',exportpath,'    Done.']);
if iscell(batchexportcallback) && ~isempty(batchexportcallback{1})
    batchexportcallback{2}=[batchexportcallback{2},{exportpath}];
end
%% Callback
callbackfun = callback{1};
callbackarg = callback{2};
callbackresult = [];
if ~isempty(callbackfun) && ~strcmp(sourceformat,'Unknown')
    disp(['Applying Callback:    ',callbackfun,'    ...']);
    callbackresult = NeuroAnalysis.Base.EvalFun(callbackfun,[{dataset},callbackarg]);
    disp(['Applying Callback:    ',callbackfun,'    Done.']);
end
%% Prepare Metadata
result.meta =[];
if ~any(strcmp(sourceformat,{'Unknown','Phy'}))
    meta = NeuroAnalysis.Base.EvalFun(...
        ['NeuroAnalysis.',sourceformat,'.PrepareMetadata'], ...
        {dataset,callbackresult});
    if (~isempty(meta))
        result.meta = meta;
    end
end
result.status = true;
