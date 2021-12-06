function [sspike] = splitspike(spike,range)
%SPLITSPIKE split out of concat spike sorting result, in which spike times are in range[start, end).

sspike=spike;
si = find(spike.time>= range(1) & spike.time< range(2));
% set spike time 0 to the start of range
sspike.time = spike.time(si)-range(1);
sspike.template = spike.template(si);
sspike.templatescale = spike.templatescale(si);
sspike.amplitude = spike.amplitude(si);

% only info of the splited clusters are included
sspike.cluster = spike.cluster(si);
sspike.clusterid = unique(sspike.cluster,'sorted');
cidx = arrayfun(@(x)find(x==spike.clusterid),sspike.clusterid);
sspike.clustergood = spike.clustergood(cidx);

if isfield(spike,'clusterwaveform')
    sspike.clusterwaveform = spike.clusterwaveform(cidx,:);
    sspike.clusterwaveformfeature = spike.clusterwaveformfeature(cidx);
end

end

