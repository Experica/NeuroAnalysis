function [di,dv,pd] = parsecrtdigitalinanalog(stream,pp,ht,lt)
%PARSECRTDIGITALINANALOG Get digital photodiode flips on a CRT display
%   Detailed explanation goes here

[ph,pi]=findpeaks(single(stream),'MinPeakProminence',pp);
[di,dv]=NeuroAnalysis.Base.parsedigitalinanalog(ph,length(ph),ht,lt);
di=pi(di);pd=mean(diff(pi));

end

