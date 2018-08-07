function [vlabtime] = reftime2vlabtime(reft,t0,timerdriftspeed,latency)
%REFTIME2VLABTIME Convert Reference Timing to VLab Timing
%   Detailed explanation goes here

if nargin==3
    latency=0;
end

vlabtime = (reft-t0)/(1+timerdriftspeed);
if latency~=0
    vlabtime = vlabtime - latency;
end

end

