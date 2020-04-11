function [tempcoords,spikeAmps,tempAmps,templates_maxwaveform_chidx,templates_maxwaveform,templates_maxwaveform_feature] = templatefeature(temps,winv,coords,chmaskradius,spikeTemplates,tempScalingAmps,fs)
%TEMPLATEFEATURE Summary of this function goes here
%   Detailed explanation goes here

tempsUnW = zeros(size(temps));
for t = 1:size(temps,1)
    tempsUnW(t,:,:) = squeeze(temps(t,:,:))*winv;
end

% The template amplitude on each channel is the height between trough to peak
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel
[tempAmpsUnscaled,maxch] = max(tempChanAmps,[],2);

% template center of mass from masked weighted channel positions
tempcoords = zeros(size(temps,1),size(coords,2));
for t = 1:size(temps,1)
    maxchpos = coords(maxch(t),:);
    dtomaxch = coords-maxchpos;
    spreadch = arrayfun(@(i)norm(dtomaxch(i,:)),1:size(coords,1)) <= chmaskradius;
    h = tempChanAmps(t,spreadch)';
    p = coords(spreadch,:);
    tempcoords(t,:) = sum(h.*p)/sum(h);
end

% assign all spikes the amplitude of their template multiplied by their
% scaling amplitudes (templates are zero-indexed)
spikeAmps = tempAmpsUnscaled(spikeTemplates+1).*tempScalingAmps;

% take the average of all spike amps to get actual template amps (since
% tempScalingAmps are equal mean for all templates)
ta = clusterAverage(spikeTemplates+1, spikeAmps);
tids = unique(spikeTemplates);
tempAmps(tids+1) = ta; % because ta only has entries for templates that had at least one spike
tempAmps = tempAmps'; % for consistency, make first dimension template number

% Get channel with largest amplitude, take that as the waveform
[~,max_ch] = max(max(abs(temps),[],2),[],3);
templates_maxwaveform = nan(size(temps,1),size(temps,2));
templates_maxwaveform_feature = struct([]);
for t = 1:size(temps,1)
    templates_maxwaveform(t,:) = temps(t,:,max_ch(t));
    templates_maxwaveform_feature = [templates_maxwaveform_feature;NeuroAnalysis.Base.spikefeature(templates_maxwaveform(t,:),fs)];
end

templates_maxwaveform_chidx = max_ch;

end

