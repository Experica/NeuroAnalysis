function [reft] = time2ref(time,t0,timerdriftspeed,latency)
%TIME2REF Convert Time to Reference Time
%   Detailed explanation goes here

if nargin==3
    latency=0;
end

reft = time*(1+timerdriftspeed)+t0+latency;

end

