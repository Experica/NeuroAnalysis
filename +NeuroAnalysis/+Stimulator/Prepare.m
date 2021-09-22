function [ stimulatordataset ] = Prepare( filepath,varargin )
%PREPARE Read and Prepare Stimulator experimental data
%   Detailed explanation goes here

p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addOptional(p,'dataset',struct([]));
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
dataset = p.Results.dataset;
%% Check file
stimulatordataset = [];
if(exist(filepath,'file')~=2)
    warning(['File do not exist: ',filepath]);
    return;
end
d = dir(filepath);
%% Read data
% disp(['Reading Stimulator File:    ',filepath,'    ...']);
% ex = dataset.ex;
% disp(['Reading StimulatorFile:    ',filepath,'    Done.']);
%% Prepare data
if ~isempty(dataset)
    if strcmp(dataset.DAQformat,'Scanbox') && (strcmp(dataset.source(1:3),'AF9') ||strcmp(dataset.source(1:3),'AF7') || strcmp(dataset.source(1:3),'AF4') || strcmp(dataset.source(1:3),'AF3') || strcmp(dataset.source(1:3),'AF2') || strcmp(dataset.source(1:3),'AF1') || strcmp(dataset.source(1:3),'AE7')|| strcmp(dataset.source(1:3),'AE6')|| strcmp(dataset.source(1:3),'AE5'))
        dataset.Version=0;
    else
        dataset.Version=1;
    end
    disp(['Preparing Stimulator Data(v',num2str(dataset.Version),'):    ',filepath,'    ...']);
    stimulatordataset = NeuroAnalysis.Base.EvalFun(['NeuroAnalysis.Stimulator.prepare',num2str(dataset.Version)],{filepath,varargin});
    
    stimulatordataset.sourceformat = 'Stimulator';
    stimulatordataset.DAQformat = 'Scanbox';
    stimulatordataset.datenum = d.datenum;
    stimulatordataset.date = d.date;  
%     stimulatordataset.source = filepath;
%     stimulatordataset.sourceformat = 'Stimulator';
    disp(['Preparing Stimulator Data(v',num2str(dataset.Version),'):    ',filepath,'    Done.']);
end

end