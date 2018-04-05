function [ visstimdataset ] = Prepare( filepath, varargin )
%PREPARE Read and Prepare VisStim data
%   Detailed explanation goes here

import NeuroAnalysis.VisStim.*

p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addOptional(p,'dataset',struct([]));
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
dataset = p.Results.dataset;
%% check file
visstimdataset = struct([]);
if ~exist(filepath, 'file')
    error(['Can not open file: ',filepath]);
end
d = dir(filepath);
%% Read data
disp(['Reading VisStim File:    ',filepath,'    ...']);
ex = struct;
ex.raw = load(filepath, '-regexp', '^(?!handles)^(?!hObject)\w');
ex.source = filepath;
ex.sourceformat = 'VisStim';
disp('Reading VisStim File:    Done.');
%% Prepare data
disp(['Preparing VisStim File:    ',filepath,'    ...']);
if isempty(ex.raw)
    warning('Preparing VisStim File:    Empty datafile');
    return;
end
ex = parseex( ex );
% Determine ripple start time
ex.t0=0;
if ~isempty(dataset)
    if isfield(dataset,'digital')
        startdchidx = find(arrayfun(@(x)x.channel==visstimconfig.StartDCh,dataset.digital));
        if ~isempty(startdchidx)
            ex.t0=dataset.digital(startdchidx).time(1);
        end
    end
end
% Fix stim times according to new start and the config file
ex = NeuroAnalysis.VisStim.adjustStimTimes(ex,...
    dataset,visstimconfig);
% Standardize experimental parameters
ex = NeuroAnalysis.Base.StandardizeEx(ex);
% Organize into dataset
visstimdataset = ex;
visstimdataset.date = d.datenum;

disp('Preparing VisStim File:    Done.');

end
