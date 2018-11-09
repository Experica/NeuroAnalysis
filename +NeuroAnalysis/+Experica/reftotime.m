function [t] = reftotime(reft,t0,timerdriftspeed,latency)
%REFTOTIME Convert Reference Time to Original Time
%   Detailed explanation goes here

if nargin==3
    latency=0;
end

t = (reft-t0-latency)/(1+timerdriftspeed);

end

