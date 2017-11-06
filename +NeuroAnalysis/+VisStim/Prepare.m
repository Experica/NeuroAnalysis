function [ dataset ] = Prepare( filepath, varargin )
%PREPARE Read and Prepare VisStim data
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'filepath');
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
%% check file
dataset = struct([]);
if ~exist(filepath, 'file') || isempty(whos('-file', filepath))
    error(['Can not open file: ',filepath]);
end
dataset = struct;
%% Read data
disp(['Reading VisStim File:    ',filepath,'    ...']);
fileVars = who('-file', filepath);
if ismember('Params',fileVars)
    load(filepath, 'Params')
    if isfield(Params, 'Data') && ~isempty(Params.Data)
        dataset.ex = Params;
    end
end
if isempty(fieldnames(dataset)) % convert old data
    [datapath, filename, ~] = fileparts(filepath);
    dataset.ex = convertLogFile( datapath, filename );
end
if ~isempty(fieldnames(dataset))
    dataset.source = filepath;
    dataset.sourceformat = 'VisStim';
end
disp('Reading VisStim File:    Done.');
%% Prepare data

end
