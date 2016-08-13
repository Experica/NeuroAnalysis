function [ result ] = Export(dataset, exportpath,varargin )
%EXPORT Export prepared dataset in Matlab MAT file
%   Detailed explanation goes here

%% Parse arguments
p = inputParser;
addRequired(p,'dataset');
addRequired(p,'exportpath');
addOptional(p,'sourceformat','Ripple',@(x)isa(x,'char'));
addOptional(p,'isparallel',true,@(x)isa(x,'logical'));
parse(p,dataset,exportpath,varargin{:});
dataset = p.Results.dataset;
exportpath = p.Results.exportpath;
sourceformat = p.Results.sourceformat;
isparallel = p.Results.isparallel;
%% Batch export
if isa(dataset,'cell')
    funlist=repelem({'NeuroAnalysis.IO.Export'},length(dataset));
    vararginlist = arrayfun(@(i)[i,{exportpath},varargin],dataset,'UniformOutput',false);
    result = NeuroAnalysis.Base.ApplyFunctions(funlist,vararginlist,isparallel);
    return;
end
%% Get dataset and export path
result.status = false;
result.source = '';
if isa(dataset,'char')
    eval(['dataset = NeuroAnalysis.',sourceformat,'.Prepare(dataset);']);
end
[datapath,filename,ext] = fileparts(dataset.source);
exportpath = fullfile(exportpath,[filename, '.mat']);
%% Save dataset
disp(['Exporting Dataset:    ',exportpath,'    ...']);
save(exportpath,'dataset','-v7.3');
result.status = true;
result.source = filename;
disp('Exporting Dataset:    Done.');
end

