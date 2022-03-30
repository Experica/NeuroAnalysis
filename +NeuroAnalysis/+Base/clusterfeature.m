function [cluwaveforms,clucoords,clumaxwaveform,cluwaveformfeature] = clusterfeature(mmapbinfile,spiketime,spikeclu,cluid,chmap,winv,coords,chmaskradius,fs)
%CLUSTERFEATURE Summary of this function goes here
%   Detailed explanation goes here

% maximum number of waveform to average
nw=1000;
% duration before and after spike time
halfwaveduration=0.00134; % 1.34ms
halfwaven = floor(halfwaveduration*fs);
spikerange = -halfwaven:halfwaven;
ns = length(spikerange);
nch = length(chmap);
nc = length(cluid);

%% extract cluster waveform
cluwaveforms = zeros(nc,nch,ns);
sidx=[];
sci=[];
for i=1:nc
    csidx = find(spikeclu==cluid(i));
    cnw = min(nw,length(csidx)-1); % exclude the last spike
    sidx = [sidx;csidx(randperm(length(csidx)-1,cnw))];
    sci = [sci;i*ones(cnw,1)];
end

st = spiketime(sidx);
[~,isort]=sort(st);
sr = arrayfun(@(i)i+spikerange,st(isort),'uniformoutput',false);
sci=sci(isort);
waveforms = reshape(mmapbinfile.Data.ap(chmap,[sr{:}]),nch,ns,[]);

uci=unique(sci);
for i=1:length(uci)
    ci = uci(i);
    cluwaveforms(ci,:,:) = mean(waveforms(:,:,sci==ci),3,'double');
end
cluwaveforms = permute(cluwaveforms,[1,3,2]); % nClusters x nTimePoints x nChannels

%% cluster feature
cluwaveformsUnW = zeros(size(cluwaveforms));
for c = 1:size(cluwaveforms,1)
    cluwaveformsUnW(c,:,:) = squeeze(cluwaveforms(c,:,:))*winv;
end

% The waveform amplitude on each channel is the height between trough and peak
cluChanAmps = squeeze(max(cluwaveformsUnW,[],2))-squeeze(min(cluwaveformsUnW,[],2));

% The cluster amplitude is the amplitude of its largest channel
[maxamp,maxch] = max(cluChanAmps,[],2);

% soma position as the cluster amplitude weighted channel positions(center of mass) within a radius centered on max cluster amplitude channel
clucoords = zeros(size(cluwaveforms,1),size(coords,2));
for c = 1:size(cluwaveforms,1)
    maxchpos = coords(maxch(c),:);
    dtomaxch = coords-maxchpos;
    spreadch = arrayfun(@(i)norm(dtomaxch(i,:)),1:size(coords,1)) <= chmaskradius;
    h = cluChanAmps(c,spreadch)';
    p = coords(spreadch,:);
    clucoords(c,:) = sum(h.*p)/sum(h);
end

clumaxwaveform = zeros(size(cluwaveforms,1),size(cluwaveforms,2));
cluwaveformfeature = cell(size(cluwaveforms,1),1);
for c = 1:size(cluwaveforms,1)
    clumaxwaveform(c,:) = cluwaveformsUnW(c,:,maxch(c)); % waveform from largest amplitude
    cluwaveformfeature{c} = NeuroAnalysis.Base.spikefeature2(squeeze(cluwaveformsUnW(c,:,:)),cluChanAmps(c,:),maxch(c),coords,fs);
end
cluwaveformfeature = [cluwaveformfeature{:}];

end

