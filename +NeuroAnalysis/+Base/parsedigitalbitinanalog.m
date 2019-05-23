function [dchs] = parsedigitalbitinanalog(stream,n,nbitpersample)
%PARSEDIGITALBITINANALOG Get digital flips in stream sample bits
%   Detailed explanation goes here

p=stream(1);pbit=de2bi(p,nbitpersample);ts=cell(nbitpersample,1);vs=cell(nbitpersample,1);
for i=2:n
    c=stream(i);
    if c==p
        continue;
    end
    cbit=de2bi(c,nbitpersample);
    for j=1:nbitpersample
        if pbit(j)~=cbit(j)
            ts{j}=[ts{j},i];
            vs{j}=[vs{j},cbit(j)];
        end
    end
    p=c;
    pbit=cbit;
end
k=1;
for i=1:nbitpersample
    t=ts{i};
    if ~isempty(t)
        v=vs{i};
        dchs(k).channel=i;
        dchs(k).time=t;
        dchs(k).value=v;
        k=k+1;
    end
end

end

