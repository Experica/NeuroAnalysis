function [reft] = toreftime(time,t0,timerdriftspeed,latency)
%TOREFTIME Convert to Reference Timing
%   Detailed explanation goes here

if nargin==3
    latency=0;
end

reft = time*(1+timerdriftspeed)+t0+latency;

end