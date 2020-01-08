function [t] = sample2time(sampleindex,fs,secondperunit,zerotimesampleindex)
%SAMPLE2TIME Summary of this function goes here
%   Detailed explanation goes here
% This is the interval time between two adjacent samples. Peichao

if nargin==2
    secondperunit=1;
    zerotimesampleindex=1;
elseif nargin==3
    zerotimesampleindex=1;
end

t=double(sampleindex-zerotimesampleindex)/fs/secondperunit;

end

