function [ex] = prepare0(ex,dataset)
%PREPARE0 Prepare version 0 of Experica.Command data format
%   Detailed explanation goes here

ex.t0=0;
if ~isempty(dataset) && isfield(dataset,'digital')
    startsyncchidx = find(arrayfun(@(x)x.channel==3,dataset.digital));
    if ~isempty(startsyncchidx)
        ex.t0=dataset.digital(startsyncchidx).time;
    end
end

ex = NeuroAnalysis.Experica.parseex(ex);
[condon,condoff] =NeuroAnalysis.Experica.parsecondonoff(ex,dataset,1,2,20);
if ~isempty(condon)
    ex.CondTest.CondOn = condon;
end
if ~isempty(condoff)
    ex.CondTest.CondOff = condoff;
end

ex = NeuroAnalysis.Base.StandardizeEx(ex);
end