function [ exdataset ] = Prepare( filepath, varargin )
%PREPARE Read and Prepare Experica.Commamd data
%   Detailed explanation goes here

p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addOptional(p,'dataset',struct([]));
addParameter(p,'exportdir','')
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
dataset = p.Results.dataset;
exportdir = p.Results.exportdir;
%% Check file
exdataset = [];
if(exist(filepath,'file')~=2)
    warning(['File do not exist: ',filepath]);
    return;
end
d = dir(filepath);
%% Read data
disp(['Reading Experica File:    ',filepath,'    ...']);
ex = yaml.ReadYaml(filepath);
disp(['Reading Experica File:    ',filepath,'    Done.']);
%% Prepare data
if ~isempty(ex)
    if ~isfield(ex,'Version')
        ex.Version=0;
    end
    disp(['Preparing Experica Data(v',num2str(ex.Version),'):    ',filepath,'    ...']);
    exdataset.ex = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.Experica.prepare',num2str(ex.Version)],{ex,dataset});
    
    exdataset.ex.source = filepath;
    exdataset.ex.sourceformat = 'Experica';
    exdataset.ex.date = d.datenum;
    
    exdataset.source = filepath;
    exdataset.sourceformat = 'Experica';
    disp(['Preparing Experica Data(v',num2str(ex.Version),'):    ',filepath,'    Done.']);
end
end
