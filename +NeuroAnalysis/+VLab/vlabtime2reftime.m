function [reft] = vlabtime2reftime(ex,vlabtime,isadddisplaylatency)
%VLABTIME2REFTIME Convert VLab Timing to Reference Timing
%   Detailed explanation goes here

if nargin==2
    isadddisplaylatency=false;
end

reft = vlabtime*(1+ex.TimerDriftSpeed)+ex.t0;
if isadddisplaylatency
    reft = reft + ex.DisplayLatency;
end

end

