function [on,off] = parsecondonoff(ex,dataset,MarkDCh,msr)
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
        dintime = dataset.digital(markdchidx).time;
        dinvalue = dataset.digital(markdchidx).data;
        if vlabconfig.NMarkPerCond*nct == length(dintime) && all(diff(dinvalue))
            isdinmarkerror = false;
        end
    end
else
    isdinmark=false;
end

iscondstate = isfield(ex.CondTest,'CONDSTATE');
for i=1:nct
    if isdinmark
        if ~isdinmarkerror
            on(i) = dintime((i-1)*vlabconfig.NMarkPerCond+1);
            off(i) = dintime((i-1)*vlabconfig.NMarkPerCond+2);
        elseif iscondstate
            tss = vlabtimetodatatime(ex,findstatetime(ex.CondTest.CONDSTATE{i},'COND'),true);
            tes = vlabtimetodatatime(ex,findstatetime(ex.CondTest.CONDSTATE{i},'SUFICI'),true);
            [isfound,ts,te]=trysearchmarktime(tss,tes,dintime,dinvalue,msr);
            if isfound
                on=[on,ts];
                off=[off,te];
            else
                on=[on,tss];
                off=[off,tes];
            end
        end
    elseif iscondstate
        on=[on,vlabtimetodatatime(ex,findstatetime(ex.CondTest.CONDSTATE{i},'COND'),true)];
        off=[off,vlabtimetodatatime(ex,findstatetime(ex.CondTest.CONDSTATE{i},'SUFICI'),true)];
    end
end

% None-ICI Mark Mode
if ex.PreICI==0 && ex.SufICI==0 && ~isempty(on) && ~isempty(off)
    for i=1:nct-1
        nextontime = on(i+1);
        currentontime=on(i);
        if (nextontime - currentontime) > (ex.CondDur+2*msr)
            off(i) = currentontime + ex.CondDur;
        else
            off(i)=nextontime;
        end
    end
    off(end)=on(end)+ex.CondDur;
end

end

