function [ vlabdataset ] = Prepare( filepath, varargin )
%PREPARE Read and Prepare VLab data
%   Detailed explanation goes here

import NeuroAnalysis.VLab.*
p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addOptional(p,'dataset',struct([]));
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
dataset = p.Results.dataset;
%% check file
vlabdataset = struct([]);
[hfile] = fopen(filepath,'r');
if hfile == -1
    warning(['Can not open file: ',filepath]);
    return;
end
vlabdataset = struct;
%% Read data
disp(['Reading VLab File:    ',filepath,'    ...']);
vlabdataset.ex = yaml.ReadYaml(filepath);
disp('Reading VLab File:    Done.');
%% Prepare data
disp(['Preparing VLab File:    ',filepath,'    ...']);
if ~isempty(vlabdataset)
    vlabdataset.ex.t0=0;
    if ~isempty(dataset)
        if isfield(dataset,'digital')
            startdchidx = find(arrayfun(@(x)x.channel==vlabconfig.StartDCh,dataset.digital));
            if ~isempty(startdchidx)
                vlabdataset.ex.t0=dataset.digital(startdchidx).time;
            end
        end
    end
    [condon,condoff] =NeuroAnalysis.VLab.parsecondonoff(vlabdataset.ex,dataset,vlabconfig.MarkDCh,vlabconfig.MarkSearchRadius);
    if ~isempty(condon)
        vlabdataset.ex.CondTest.CondOn = condon;
    end
    if ~isempty(condoff)
        vlabdataset.ex.CondTest.CondOff = condoff;
    end
    vlabdataset.ex = NeuroAnalysis.VLab.parseex(vlabdataset.ex);
    
    vlabdataset.source = filepath;
    vlabdataset.sourceformat = 'VLab';
end
disp('Preparing VLab File:    Done.');
end
