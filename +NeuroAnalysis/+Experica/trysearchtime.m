function [st] = trysearchtime(t,ints,sr)
%TRYSEARCHTIME Try to search the closest time in a time sequence
%   Detailed explanation goes here

tc=[];
for i =1:length(ints)
    d = ints(i)-t;
    if abs(d)<= sr
        tc=[tc,ints(i)];
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

