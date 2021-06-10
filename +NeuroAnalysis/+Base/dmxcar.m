function [data] = dmxcar(data,dmxgroup,rmdc)
%DMXCAR Demuxed CAR
%   Detailed explanation goes here

if nargin==2
    rmdc=false;
end

for g = 1:length(dmxgroup)
    if rmdc
        data(dmxgroup{g},:) = data(dmxgroup{g},:) - mean(data(dmxgroup{g},:),2,'native');
    end
    data(dmxgroup{g},:) = data(dmxgroup{g},:) - median(data(dmxgroup{g},:),1);
end

end