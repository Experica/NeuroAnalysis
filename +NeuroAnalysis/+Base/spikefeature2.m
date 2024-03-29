function [swf] = spikefeature2(waveforms,amps,maxch,chcoords,fs)
%SPIKEFEATURE2 2D Spike Waveforms Feature
%   Detailed explanation goes here

% percentage of max amplitude to threshold spatial spread
spreadampthr = 0.15;
% max possible spread radius(um)
maxspreadradius = 350;
% left/right distance(um) to max channel, within which are considered for up/down propagation
updownhalfwidth = 20;

maxamp = amps(maxch);
chtomaxch = chcoords-chcoords(maxch,:);
chmask = (amps >= maxamp*spreadampthr) & (arrayfun(@(i)norm(chtomaxch(i,:)),1:size(chtomaxch,1)) <= maxspreadradius);
maskedchs = find(chmask);
maxmaskedch = find(maskedchs==maxch);

chwavefeature = cell(length(maskedchs),1);
for i=1:length(maskedchs)
    chwavefeature{i} = NeuroAnalysis.Base.spikefeature(waveforms(:,maskedchs(i)),fs);
end
chwavefeature = [chwavefeature{:}];


dtomaxch = chtomaxch(chmask,:);
upchidx = dtomaxch(:,2)>=0 & abs(dtomaxch(:,1)) < updownhalfwidth;
downchidx = dtomaxch(:,2)<=0 & abs(dtomaxch(:,1)) < updownhalfwidth;
% 2D spatial spread centered at max amplitude channel
upspread = max(dtomaxch(:,2));
downspread = min(dtomaxch(:,2));
leftspread = min(dtomaxch(:,1));
rightspread = max(dtomaxch(:,1));

% inverse of upward propagation speed
y = [chwavefeature(upchidx).ttrough]';
X = [ones(size(y)), dtomaxch(upchidx,2)];
c = X\y;
uppvinv = c(2);

% inverse of downward propagation speed
y = [chwavefeature(downchidx).ttrough]';
X = [ones(size(y)), dtomaxch(downchidx,2)];
c = X\y;
downpvinv = c(2);

swf = chwavefeature(maxmaskedch); % 1D feature chosen from max amplitude channel
swf.upspread = upspread;
swf.downspread = downspread;
swf.leftspread = leftspread;
swf.rightspread = rightspread;
swf.uppvinv = uppvinv;
swf.downpvinv = downpvinv;

end

