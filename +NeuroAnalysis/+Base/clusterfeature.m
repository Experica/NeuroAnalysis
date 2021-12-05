function [cluwaveform,cluwaveformfeature] = clusterfeature(mmapbinfile,spiketime,spikeclu,cluid,chmap,chcoords,fs)
%CLUSTERFEATURE Summary of this function goes here
%   Detailed explanation goes here

% percentage of max amplitude to threshold spatial spread
spreadampthr = 0.15;
% maximum number of waveform to average
nw=2000;
% duration before and after trough point
halfwaveduration=0.00134; % 1.34ms
halfwavens = floor(halfwaveduration*fs);
spikerange = -halfwavens:halfwavens;
ns = length(spikerange);
nch = length(chmap);
cluwaveform = zeros(length(cluid),ns);
cluwaveformfeature = cell(length(cluid),1);

for i=1:length(cluid)
    cluidx = find(spikeclu==cluid(i));
    cnw = min(nw,length(cluidx)-1); % exclude the last spike
    sidx = sort(cluidx(randperm(length(cluidx)-1,cnw)));
    sr = arrayfun(@(i)spiketime(i)+spikerange,sidx,'uniformoutput',false);
    cluwaveforms = mmapbinfile.Data.ap(chmap,[sr{:}]);
    cluwaveforms = mean(reshape(cluwaveforms,nch,ns,[]),3,'double')'; % nsample x nch
    
    amp = max(cluwaveforms,[],1)-min(cluwaveforms,[],1);
    [maxamp,maxch] = max(amp);
    spreadchmask = amp >= maxamp*spreadampthr;
    
    cluwaveform(i,:) = cluwaveforms(:,maxch);
    cluwaveformfeature{i} = NeuroAnalysis.Base.spikefeature2(cluwaveforms,spreadchmask,maxch,chcoords,fs);
end
cluwaveformfeature = [cluwaveformfeature{:}];

end

