function [on,off] = parsecondonoff(ex,dataset,CondDCh,MarkDCh,lsr)
%PARSECONDONOFF Try to get the most accurate timing of Condition On/Off
%   Detailed explanation goes here

import NeuroAnalysis.Experica.*
nct = length(ex.CondTest.CondIndex);
on = [];off=[];
if ~isempty(dataset) && isfield(dataset,'digital')
    markdchidx = find(arrayfun(@(x)x.channel==MarkDCh,dataset.digital));
    isdinmark = ~isempty(markdchidx);
    isdinmarkerror=true;
    if isdinmark
        dinmarktime = dataset.digital(markdchidx).time;
        dinmarkvalue = dataset.digital(markdchidx).data;
        if 2*nct == length(dinmarktime) && all(diff(dinmarkvalue))
            isdinmarkerror = false;
        end
    end
    
    conddchidx = find(arrayfun(@(x)x.channel==CondDCh,dataset.digital));
    isdincond = ~isempty(conddchidx);
    isdinconderror=true;
    if isdincond
        dincondtime = dataset.digital(conddchidx).time;
        dincondvalue = dataset.digital(conddchidx).data;
        if 2*nct == length(dincondtime) && all(diff(dincondvalue))
            isdinconderror = false;
        end
    end
else
    isdinmark=false;
    isdincond=false;
end

    function [isfound,ts,te] = trysearchmarktime(tss,tes,dintime,dinvalue,sr)
        %TRYSEARCHMARKTIME Try to search mark time
        
        msv = dinvalue(1);
        tssidx=[];tesidx=[];
        for j=1:length(dinvalue)
            dts = dintime(j)-tss;
            dte = dintime(j)-tes;
            if (dinvalue(j) == msv)
                if(abs(dts)<=sr)
                    tssidx=[tssidx,j];
                end
            else
                if(abs(dte)<=sr)
                    tesidx=[tesidx,j];
                end
            end
            if dts>sr&&dte>sr
                break;
            end
        end
        if length(tssidx)==1 && length(tesidx)==1
            ts=dintime(tssidx(1));
            te=dintime(tesidx(1));
            isfound=true;
        else
            ts=0;te=0;isfound=false;
        end
        
    end

for i=1:nct
    if isdinmark
        if ~isdinmarkerror
            on(i) = dinmarktime((i-1)*2+1);
            off(i) = dinmarktime((i-1)*2+2);
        else
            if isdincond && ~isdinconderror % use digital in condition times if possible
                tss = dincondtime((i-1)*2+1) + ex.Latency;
                tes = dincondtime((i-1)*2+2) + ex.Latency;
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
        on(i) = dincondtime((i-1)*2+1) + ex.Latency;
        off(i) = dincondtime((i-1)*2+2) + ex.Latency;
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