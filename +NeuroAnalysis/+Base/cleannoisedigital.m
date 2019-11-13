function [cdt,cdv] = cleannoisedigital(dt,dv,minlowdur,minhighdur)
%CLEANNOISEDIGITAL Summary of this function goes here
%   Detailed explanation goes here

cdt=dt(1);cdv=dv(1);
hi=find(dv==1);
for i=2:length(hi)
    chi = hi(i);
    cli = hi(i)-1;
    clt = dt(cli);
    cht = dt(chi);
    if ((cht-clt) >= minlowdur) && ((clt-cdt(end)) >= minhighdur)
        cdt = [cdt,clt,cht];
        cdv = [cdv,dv(cli),dv(chi)];
    end
end
if dv(end)==0
    cdt=[cdt,dt(end)];
    cdv=[cdv,0];
end

end

