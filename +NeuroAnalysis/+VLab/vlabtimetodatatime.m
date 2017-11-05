function [datatime] = vlabtimetodatatime(ex,vlabtime,isaddlatency)
%VLABTIMETODATATIME Convert VLab Timing to Data Timing
%   Detailed explanation goes here

if nargin==2
    isaddlatency=false;
end

datatime = vlabtime*(1+ex.TimerDriftSpeed)+ex.t0;
if isaddlatency
    datatime = datatime + ex.Latency;
end

end

