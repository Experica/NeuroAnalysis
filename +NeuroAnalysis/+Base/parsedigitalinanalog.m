function [di,dv] = parsedigitalinanalog(stream,n,ht,lt)
%PARSEDIGITALINANALOG Get digital rising and falling edges in a analog stream
%   Detailed explanation goes here

di=[];dv=[];cs=stream(1)>ht;
for i=2:n
    t=stream(i);
    if t > ht && ~cs
        cs=true;
        di=[di,i];
        dv=[dv,cs];
    elseif t < lt && cs
        cs=false;
        di=[di,i];
        dv=[dv,cs];
    end
end

end

