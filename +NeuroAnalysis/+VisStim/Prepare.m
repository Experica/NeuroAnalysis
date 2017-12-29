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
visstimdataset = struct;
%% Read data
disp(['Reading VisStim File:    ',filepath,'    ...']);
visstimdataset.ex = struct;
visstimdataset.ex.DataPath = filepath;
visstimdataset.ex.raw = load(filepath, '-regexp', '^(?!handles)^(?!hObject)\w');
disp('Reading VisStim File:    Done.');
%% Prepare data
disp(['Preparing VisStim File:    ',filepath,'    ...']);
if ~isempty(visstimdataset.ex.raw)
    visstimdataset.ex = parseex( visstimdataset.ex );
    visstimdataset.ex.t0=0;
    if ~isempty(dataset)
        if isfield(dataset,'digital')
            startdchidx = find(arrayfun(@(x)x.channel==visstimconfig.StartDCh,dataset.digital));
            if ~isempty(startdchidx)
                visstimdataset.ex.t0=dataset.digital(startdchidx).time(1);
            end
        end
    end
    visstimdataset.ex = NeuroAnalysis.VisStim.adjustStimTimes(visstimdataset.ex,...
        dataset,visstimconfig);
    
    visstimdataset.source = filepath;
    visstimdataset.sourceformat = 'VisStim';
end
disp('Preparing VisStim File:    Done.');

end
