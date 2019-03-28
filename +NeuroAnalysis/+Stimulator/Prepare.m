function [ stimulatordataset ] = Prepare( filepath,varargin )
%PREPARE Read and Prepare Stimulator data
%   Detailed explanation goes here

import NeuroAnalysis.Stimulator.*

p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addOptional(p,'dataset',struct([]));
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
dataset = p.Results.dataset;
%% check file
stimulatordataset = struct([]);
if ~exist(filepath, 'file')
    error(['Can not open file: ',filepath]);
end
d = dir(filepath);
%% Read data
disp(['Reading Stimulator File:    ',filepath,'    ...']);
ex = struct;
ex.raw = load(filepath,'-mat');
ex.raw = ex.raw.Analyzer;
ex.source = filepath;
ex.sourceformat = 'Stimulator';
disp('Reading Stimulator File:    Done.');
%% Prepare data
disp(['Preparing Stimulator File:    ',filepath,'    ...']);
if isempty(ex.raw)
    warning('Preparing Stimulator File:    Empty datafile');
    return;
end
ex = parseex( ex );
% Standardize experimental parameters
ex = NeuroAnalysis.Base.StandardizeEx(ex);
% Organize into dataset
stimulatordataset = ex;
stimulatordataset.date = d.datenum;

disp('Preparing Stimulator File:    Done.');
end

