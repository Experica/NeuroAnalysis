function [on,off] = parsecondonoff(ex,dataset,CondDCh,MarkDCh,lsr)
%PARSECONDONOFF Try to get the most accurate timing of Condition On/Off
%   Detailed explanation goes here

import NeuroAnalysis.VLab.*
nct = length(ex.CondTest.CondIndex);
on = [];off=[];
if ~isempty(dataset) && isfield(dataset,'digital')
    markdchidx = find(arrayfun(@(x)x.channel==MarkDCh,dataset.digital));
    isdinmark = ~isempty(markdchidx);
    isdinmarkerror=true;
    if isdinmark
        dinmarktime = dataset.digital(markdchidx).time;
        dinmarkvalue = dataset.digital(markdchidx).data;
        if vlabconfig.NMarkPerCond*nct == length(dinmarktime) && all(diff(dinmarkvalue))
            isdinmarkerror = false;
        end
    end
    
    conddchidx = find(arrayfun(@(x)x.channel==CondDCh,dataset.digital));
    isdincond = ~isempty(conddchidx);
    isdinconderror=true;
    if isdincond
        dincondtime = dataset.digital(conddchidx).time;
        dincondvalue = dataset.digital(conddchidx).data;
        if vlabconfig.NMarkPerCond*nct == length(dincondtime) && all(diff(dincondvalue))
            isdinconderror = false;
        end
    end
else
    isdinmark=false;
    isdincond=false;
end

for i=1:nct
    if isdinmark
        if ~isdinmarkerror
            on(i) = dinmarktime((i-1)*vlabconfig.NMarkPerCond+1);
            off(i) = dinmarktime((i-1)*vlabconfig.NMarkPerCond+2);
        else
            if isdincond && ~isdinconderror % use digital in condition times if possible
                tss = dincondtime((i-1)*vlabconfig.NMarkPerCond+1) + ex.Latency;
                tes = dincondtime((i-1)*vlabconfig.NMarkPerCond+2) + ex.Latency;
            else
                tss = ex.CondTest.CondOnTime(i);
                tes = ex.CondTest.SufICIOnTime(i);
            end
            [isfound,ts,te]=trysearchmarktime(tss,tes,dinmarktime,dinmarkvalue,lsr);
            if isfound
                on(i)=ts;
                off(i)=te;
            else
                on(i)=tss;
                off(i)=tes;
            end
        end
    elseif isdincond && ~isdinconderror
        on(i) = dincondtime((i-1)*vlabconfig.NMarkPerCond+1) + ex.Latency;
        off(i) = dincondtime((i-1)*vlabconfig.NMarkPerCond+2) + ex.Latency;
    else
        on(i) = ex.CondTest.CondOnTime(i);
        off(i) = ex.CondTest.SufICIOnTime(i);
    end
end

% None-ICI Mark Mode
if ex.PreICI==0 && ex.SufICI==0 && ~isempty(on) && ~isempty(off)
    for i=1:nct-1
        currentontime=on(i);
        nextontime = on(i+1);
        if (nextontime - currentontime) > (ex.CondDur+2*lsr)
            off(i) = currentontime + ex.CondDur;
        else
            off(i)=nextontime;
        end
    end
    off(end)=on(end)+ex.CondDur;
end

end

