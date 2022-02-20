function [spike] = extractrez3(rez,secondperunit)
%EXTRACTREZ3 extract spike directly from kilosort3 result
%   Detailed explanation goes here

if nargin==1
    secondperunit=1;
end

chmaskradius = 55; % radius(um) within which the templates height are used to estimate position

W = gather(single(rez.W));
U = gather(single(rez.U));
nt0 = size(W,1);
Nfilt = size(W,2);
Nch = rez.ops.Nchan;
coords = [rez.xcoords(:) rez.ycoords(:)];
chanMap = rez.ops.chanMap(:);
spikeTimes = rez.st3(:,1);
spikeTemplates = rez.st3(:,2);
tempScalingAmps = rez.st3(:,3);

% each template in spatial(electrodes) and temporal(waveform) dimention
temps = zeros(Nch, nt0, Nfilt, 'single');
for iNN = 1:Nfilt
    temps(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
end
temps = permute(temps, [3 2 1]); % nTemplates x nSamples x nChannels

wraw = rez.Wrot/rez.ops.scaleproc;
winvraw = wraw^-1;
w = eye(size(rez.Wrot)) / rez.ops.scaleproc;
winv = w^-1;


rez.mu = gather(single(rez.mu));
rez.ops.igood = gather(rez.ops.igood);
spike.ops = rez.ops;
spike.fs = rez.ops.fs;
% use correct whiteningmatrixinv to restore unwhiten template waveform
[tempcoords,spikeAmps,tempAmps,templates_maxwaveform,templates_waveform_feature]...
    = NeuroAnalysis.Base.templatefeature(temps,winvraw,coords,chmaskradius,spikeTemplates-1,tempScalingAmps,spike.fs);
spike.templates = temps;
spike.templatesposition = tempcoords;
spike.templatesamplitude = tempAmps;
spike.templateswaveform = templates_maxwaveform;
spike.templateswaveformfeature = templates_waveform_feature;

spike.chanmap = int32(chanMap);
spike.channelposition = coords;
spike.whiteningmatrix = w;
spike.whiteningmatrixinv = winv;
spike.whiteningmatrixraw = wraw;
spike.whiteningmatrixinvraw = winvraw;

% times for each spike, t0 is the first sample in the data stream
spike.time = NeuroAnalysis.Base.sample2time(spikeTimes,spike.fs,secondperunit);
% ids of template on which each spike is extracted
spike.template = int64(spikeTemplates);
% template scaling for each spike
spike.templatescale = tempScalingAmps;
% scaled unwhiten template amplitude for each spike
spike.amplitude = spikeAmps;

% ids of cluster where each is expected to be a single cell, here
% regarding each unique template as a single cell
spike.cluster = spike.template;
% unique cluser ids(already sorted in `unique`)
spike.clusterid = unique(spike.cluster);
spike.clustergood = rez.good;


NeuroAnalysis.Visualization.plotdriftmap(spike.time,spike.amplitude,spike.templatesposition(spike.template,2));


spike.qcversion = 'kilosort3';
spike.qc = [];
spike=[];
end

