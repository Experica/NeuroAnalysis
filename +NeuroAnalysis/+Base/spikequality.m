function [spike] = spikequality(spike,isipercentile,isithreshold)
%SPIKEQUALITY Summary of this function goes here
%   Detailed explanation goes here

arguments
    spike
    % percentile of isi cumulative distribution
    isipercentile = 5
    % maximum time for isi violation
    isithreshold = 0.001/spike.secondperunit;
end

qm=struct;
for i=1:length(spike.clusterid)
    st = sort(spike.time(spike.cluster==spike.clusterid(i)));
    fr = length(st)/(st(end)-st(1));
    qm(i).fr = fr/spike.secondperunit;
    
    isi = diff(st);
    qm(i).pisi = prctile(isi,isipercentile);
    qm(i).fp = isiviolations(isi,fr,isithreshold);
end
spike.qm = qm;

    function [fp] = isiviolations(isi,fr,threshold)
        % Based on metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
        % modified by Dan Denman from cortex-lab/sortingQuality GitHub by Nick Steinmetz
        % modified by Jennifer Colonell to correctly solve the quadratic equation
        
        % fp : rate of contaminating spikes as a fraction of overall rate
        %      A perfect unit has a fp = 0
        %      A unit with some contamination has a fp < 0.5
        %      A unit with lots of contamination has a fp > 1.0
        
        nspike = length(isi)+1;
        nviolations = sum(isi<threshold);
        violationtime = 2*nspike*threshold;
        c = nviolations/(violationtime*fr);
        if c < 0.25 % valid solution of quadratic eq for fp
            fp = (1-sqrt(1-4*c))/2;
        else % no valid solution of quadratic eq, set to 1
            fp = 1;
        end
    end

end

