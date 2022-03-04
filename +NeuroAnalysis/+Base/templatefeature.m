function [tempcoords,spikeAmps,tempAmps,templates_maxwaveform,templates_waveform_feature] = templatefeature(temps,winv,coords,chmaskradius,spikeTemplates,tempScalingAmps,fs)
%TEMPLATEFEATURE Summary of this function goes here
%   Detailed explanation goes here

tempsUnW = zeros(size(temps));
for t = 1:size(temps,1)
    tempsUnW(t,:,:) = squeeze(temps(t,:,:))*winv;
end

% The template amplitude on each channel is the height between trough and peak
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel
[tempAmpsUnscaled,maxch] = max(tempChanAmps,[],2);

% soma position as the template amplitude weighted channel positions(center of mass) within a radius centered on max template amplitude channel
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
ta = NeuroAnalysis.Base.clusterAverage(spikeTemplates+1, spikeAmps);
tids = unique(spikeTemplates);
tempAmps(tids+1) = ta; % because ta only has entries for templates that had at least one spike
tempAmps = tempAmps'; % for consistency, make first dimension template number


templates_maxwaveform = zeros(size(temps,1),size(temps,2));
templates_waveform_feature = cell(size(temps,1),1);
for t = 1:size(temps,1)
    templates_maxwaveform(t,:) = tempsUnW(t,:,maxch(t)); % waveform from largest amplitude
    templates_waveform_feature{t} = NeuroAnalysis.Base.spikefeature2(squeeze(tempsUnW(t,:,:)),tempChanAmps(t,:),maxch(t),coords,fs);
end
templates_waveform_feature = [templates_waveform_feature{:}];

end

