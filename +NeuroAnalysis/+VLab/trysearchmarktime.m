function [isfound,ts,te] = trysearchmarktime(tss,tes,dintime,dinvalue,msr)
%TRYSEARCHMARKTIME Try to search mark time based on vlab time
%   Detailed explanation goes here

msv = dinvalue(1);
tssidx=[];tesidx=[];
for i=1:length(dinvalue)
    dts = dintime(i)-tss;
    dte = dintime(i)-tes;
    if(abs(dts)<=msr && dinvalue(i)==msv)
        tssidx=[tssidx,i];
    end
    if(abs(dte)<=msr && dinvalue(i) ~=msv)
        tesidx=[tesidx,i];
    end
    if dts>msr&&dte>msr
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

