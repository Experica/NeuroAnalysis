function [sspike] = splitspike(spike,range)
%SPLITSPIKE split out of concat spike sorting result, in which spike times are in range[start, end).
% only unit info about the spikes in range are included.

sspike=spike;
si = find(spike.time>= range(1) & spike.time< range(2));
% set spike time 0 to the start of range
sspike.time=spike.time(si)-range(1);
sspike.amplitude=spike.amplitude(si);
sspike.template=spike.template(si);
sspike.cluster = sspike.cluster(si);
sspike.clusterid = unique(sspike.cluster);
% only info of clusters split are included
cidx = arrayfun(@(x)find(x==spike.clusterid),sspike.clusterid);
sspike.clustertemplates = spike.clustertemplates(cidx,:,:);
sspike.clustergood = spike.clustergood(cidx);

end

