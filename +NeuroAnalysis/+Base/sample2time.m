function [t] = sample2time(sampleindex,fs,secondperunit,zerotimesampleindex)
%SAMPLE2TIME Get the times of samples, relative to zerotimesampleindex, and
% convert to time unit specified in secondperunit

if nargin==2
    secondperunit=1;
    zerotimesampleindex=1;
elseif nargin==3
    zerotimesampleindex=1;
end

t=double(sampleindex-zerotimesampleindex)/fs/secondperunit;

end

