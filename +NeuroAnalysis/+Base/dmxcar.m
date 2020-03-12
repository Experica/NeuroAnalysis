function [data] = dmxcar(data,dmxgroup)
%DMXCAR Demuxed CAR
%   Detailed explanation goes here

for g = 1:length(dmxgroup)
    data(dmxgroup{g},:) = data(dmxgroup{g},:) - mean(data(dmxgroup{g},:),2,'native');
    data(dmxgroup{g},:) = data(dmxgroup{g},:) - median(data(dmxgroup{g},:),1);
end

end