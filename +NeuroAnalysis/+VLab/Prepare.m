function [ vlabdataset ] = Prepare( filepath, varargin )
%PREPARE Read and Prepare VLab data
%   Detailed explanation goes here

p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addOptional(p,'dataset',struct([]));
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
dataset = p.Results.dataset;

%% check file
vlabdataset = [];
if(exist(filepath,'file')~=2)
    warning(['File do not exist: ',filepath]);
    return;
end
d = dir(filepath);
%% Read data
disp(['Reading VLab File:    ',filepath,'    ...']);
ex = yaml.ReadYaml(filepath);
disp(['Reading VLab File:    ',filepath,'    Done.']);
%% Prepare data
disp(['Preparing VLab Data:    ',filepath,'    ...']);
if ~isempty(ex)
    if ~isfield(ex,'Version')
        ex.Version=0;
    end
    vlabdataset.ex = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.VLab.prepare',num2str(ex.Version)],{ex,dataset});
    
    vlabdataset.ex.source = filepath;
    vlabdataset.ex.sourceformat = 'VLab';
    vlabdataset.ex.date = d.datenum;
    
    vlabdataset.source = filepath;
    vlabdataset.sourceformat = 'VLab';
end
disp(['Preparing VLab Data:    ',filepath,'    Done.']);
end
