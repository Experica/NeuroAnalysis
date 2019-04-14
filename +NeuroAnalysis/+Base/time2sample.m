function [sampleindex] = time2sample(t,fs,secondperunit,zerotimesampleindex)
%TIME2SAMPLE Summary of this function goes here
%   Detailed explanation goes here

if nargin==2
    secondperunit=1;
    zerotimesampleindex=1;
elseif nargin==3
    zerotimesampleindex=1;
end

sampleindex=round(t*secondperunit*fs)+zerotimesampleindex;

end

