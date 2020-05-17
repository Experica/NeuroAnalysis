function [cdt,cdv] = cleannoisedigital(dt,dv,minlowdur,minhighdur)
%CLEANNOISEDIGITAL clean noisy digital caused by logical high very close to threshold.
% The low state won't generate any noise, so each flips around a conddur long low epoch would be the real high states.
%
% Assume init state is low and no low-low or high-high events ever occured.

cdt=dt(1);cdv=dv(1);
his=find(dv==1);
for i=2:length(his)
    hi = his(i);
    li = hi-1;
    lt = dt(li);
    ht = dt(hi);
    if ((ht-lt) >= minlowdur) && ((lt-cdt(end)) >= minhighdur)
        cdt = [cdt,lt,ht];
        cdv = [cdv,dv(li),dv(hi)];
    end
end
if dv(end)==0
    cdt=[cdt,dt(end)];
    cdv=[cdv,0];
end

end

