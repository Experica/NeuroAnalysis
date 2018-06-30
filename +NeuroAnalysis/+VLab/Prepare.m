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

import NeuroAnalysis.VLab.*
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
disp('Reading VLab File:    Done.');
%% Prepare data
disp('Preparing VLab Data:    ...');
if ~isempty(ex)
    vlabdataset.ex = ex;
    vlabdataset.ex.t0=0;
    if ~isempty(dataset)
        if isfield(dataset,'digital')
            startdchidx = find(arrayfun(@(x)x.channel==vlabconfig.StartDCh,dataset.digital));
            if ~isempty(startdchidx)
                vlabdataset.ex.t0=dataset.digital(startdchidx).time;
            end
        end
    end
    vlabdataset.ex = NeuroAnalysis.VLab.parseex(vlabdataset.ex);
    [condon,condoff] =NeuroAnalysis.VLab.parsecondonoff(vlabdataset.ex,dataset,vlabconfig.CondDCh,vlabconfig.MarkDCh,vlabconfig.MarkSearchRadius);
    if ~isempty(condon)
        vlabdataset.ex.CondTest.CondOn = condon;
    end
    if ~isempty(condoff)
        vlabdataset.ex.CondTest.CondOff = condoff;
    end
    
    vlabdataset.ex = NeuroAnalysis.Base.StandardizeEx(vlabdataset.ex);
    vlabdataset.ex.source = filepath;
    vlabdataset.ex.sourceformat = 'VLab';
    vlabdataset.ex.date = d.datenum;
    
    vlabdataset.source = filepath;
    vlabdataset.sourceformat = 'VLab';
    
end
disp('Preparing VLab Data:    Done.');
end
