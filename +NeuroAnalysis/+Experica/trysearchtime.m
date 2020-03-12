function [st] = trysearchtime(t,tseq,sr)
%TRYSEARCHTIME Try to search the closest time in a time sequence around a perdicted time within a search radius
%   Detailed explanation goes here

tc=[];
for i =1:length(tseq)
    d = tseq(i)-t;
    if abs(d)<= sr
        tc=[tc,tseq(i)];
    end
    if d> sr
        break;
    end
end

if isempty(tc) || length(tc)>1
    st=NaN;
else
    st=tc;
end

end