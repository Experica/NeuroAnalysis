function [ visstimdataset ] = Prepare( filepath, exportdir, varargin )
%PREPARE Read and Prepare VisStim data
%   Detailed explanation goes here

import NeuroAnalysis.VisStim.*

p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addRequired(p,'exportdir');
addOptional(p,'dataset',struct([]));
parse(p,filepath,exportdir,varargin{:});
filepath = p.Results.filepath;
exportdir = p.Results.exportdir;
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
    % Determine ripple start time
    visstimdataset.ex.t0=0;
    if ~isempty(dataset)
        if isfield(dataset,'digital')
            startdchidx = find(arrayfun(@(x)x.channel==visstimconfig.StartDCh,dataset.digital));
            if ~isempty(startdchidx)
                visstimdataset.ex.t0=dataset.digital(startdchidx).time(1);
            end
        end
    end
    % Fix stim times according to new start and the config file
    visstimdataset.ex = NeuroAnalysis.VisStim.adjustStimTimes(visstimdataset.ex,...
        dataset,visstimconfig);
    
    visstimdataset.source = filepath;
    visstimdataset.sourceformat = 'VisStim';
    visstimdataset.filepath = fullfile(exportdir, ...
        sprintf('%s#%d_%s_%s.mat', visstimdataset.ex.Subject_ID, ...
        visstimdataset.ex.File_ID, visstimdataset.ex.RecordSite, ...
        visstimdataset.ex.ID));
end
disp('Preparing VisStim File:    Done.');

end
