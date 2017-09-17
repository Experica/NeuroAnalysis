function [ dataset ] = Prepare( filepath, varargin )
%PREPARE Read and Prepare VLab data
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'filepath');
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
%% check file
dataset = struct([]);
[hfile] = fopen(filepath,'r');
if hfile == -1
    warning(['Can not open file: ',filepath]);
    return;
end
dataset = struct;
%% Read data
disp(['Reading VLab File:    ',filepath,'    ...']);
dataset.ex = yaml.ReadYaml(filepath);
if ~isempty(dataset)
    dataset.source = filepath;
    dataset.sourceformat = 'VLab';
end
disp('Reading VLab File:    Done.');
%% Prepare data

end
