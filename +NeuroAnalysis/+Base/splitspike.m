function [sspike] = splitspike(spike,range)
%SPLITSPIKE split concat spike sorting result into each dataset,
% with consistent unit info within dataset, but consistent cluster id across datasets

sspike=spike;
si = find(spike.time>= range(1) & spike.time< range(2));
sspike.time=spike.time(si)-range(1);
sspike.amplitude=spike.amplitude(si);
sspike.template=spike.template(si);
sspike.cluster = sspike.cluster(si);
sspike.clusterid = unique(sspike.cluster);
% match clusterid
cidx = arrayfun(@(x)find(x==spike.clusterid),sspike.clusterid);
sspike.clustertemplates = spike.clustertemplates(cidx,:,:);
sspike.clustergood = spike.clustergood(cidx);

end

