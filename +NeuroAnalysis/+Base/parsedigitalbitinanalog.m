function [dchs] = parsedigitalbitinanalog(stream,n,nbitpersample)
%PARSEDIGITALBITINANALOG Get digital flips in stream sample bits.
% Return values(0/1) and times(sample index) for each channel in which digital signal are found.

p=stream(1);pbit=bitget(p,1:nbitpersample);ts=cell(nbitpersample,1);vs=cell(nbitpersample,1);dchs=[];
for i=2:n
    c=stream(i);
    if c==p
        continue;
    end
    cbit=bitget(c,1:nbitpersample);
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
        dchs(k).data=v;
        k=k+1;
    end
end

end

