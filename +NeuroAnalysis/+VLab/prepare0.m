function [ex] = prepare0(ex,dataset)
%PREPARE0 Prepare version 0 of VLab data format
%   Detailed explanation goes here

ex.t0=0;
if ~isempty(dataset)
    if isfield(dataset,'digital')
        startsyncdchidx = find(arrayfun(@(x)x.channel==vlabconfig.StartSyncDCh,dataset.digital));
        if ~isempty(startsyncdchidx)
            ex.t0=dataset.digital(startsyncdchidx).time;
        end
    end
end

ex = NeuroAnalysis.VLab.parseex(ex);
[condon,condoff] =NeuroAnalysis.VLab.parsecondonoff(ex,dataset,vlabconfig.EventSyncDCh,vlabconfig.EventMeasureDCh,vlabconfig.MaxDisplayLatencyError);
if ~isempty(condon)
    ex.CondTest.CondOn = condon;
end
if ~isempty(condoff)
    ex.CondTest.CondOff = condoff;
end

ex = NeuroAnalysis.Base.StandardizeEx(ex);
end