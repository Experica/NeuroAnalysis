function [cluwaveform,cluwaveformfeature] = clusterfeature(mpbinfile,spiketime,spikeclu,cluid,spiketemplate,tempwaveformch,fs)
%CLUSTERFEATURE Summary of this function goes here
%   Detailed explanation goes here

nwaveform=2000;
halfwaveduration=0.00134; % 1.34ms
halfwavensample = ceil(halfwaveduration*fs);
spikerange = -halfwavensample:halfwavensample;
nss = length(spikerange);
cluwaveform = zeros(length(cluid),nss);
cluwaveformfeature = struct([]);

for i=1:length(cluid)
    cluidx = find(spikeclu==cluid(i));
    cluch = tempwaveformch(spiketemplate(cluidx(1)));
    nw = min(nwaveform,length(cluidx));
    stidx = sort(cluidx(randperm(length(cluidx),nw)));
    si = arrayfun(@(i)spiketime(i)+spikerange,stidx,'uniformoutput',false);
    cluwaveforms = mpbinfile.Data.ap(cluch,[si{:}]);
    cluwaveform(i,:) = mean(reshape(cluwaveforms,nss,[]),2);
    cluwaveformfeature = [cluwaveformfeature;NeuroAnalysis.Base.spikefeature(cluwaveform(i,:),fs)];
end

end

