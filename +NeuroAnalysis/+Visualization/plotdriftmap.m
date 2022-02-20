function plotdriftmap(spikeTimes, spikeAmps, spikeYpos)
%PLOTDRIFTMAP Summary of this function goes here
%   Detailed explanation goes here

figure;
set(gcf, 'Color', 'w')

nColorBins = 100;
ampRange = quantile(spikeAmps, [0.05 0.95]);
colorBins = linspace(ampRange(1), ampRange(2), nColorBins);

colors = gray(nColorBins); colors = colors(end:-1:1, :);
for b = 1:nColorBins-1
    theseSpikes = spikeAmps>=colorBins(b) & spikeAmps<=colorBins(b+1);
    
    plot(spikeTimes(theseSpikes), spikeYpos(theseSpikes), '.', 'Color', colors(b,:));
    hold on;
end
axis tight
box off

xlabel('time')
ylabel('spike position (um)')
title('Drift map')
drawnow
end

