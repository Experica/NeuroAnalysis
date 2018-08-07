function [reft] = vlabtime2reftime(vlabtime,t0,timerdriftspeed,latency)
%VLABTIME2REFTIME Convert VLab Timing to Reference Timing
%   Detailed explanation goes here

if nargin==3
    latency=0;
end

reft = vlabtime*(1+timerdriftspeed)+t0;
if latency~=0
    reft = reft + latency;
end

end

