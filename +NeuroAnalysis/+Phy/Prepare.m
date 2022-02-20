function [dataset] = Prepare(phydir,varargin)
%PREPARE Read Phy result in phydir and merge back to dataset
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'phydir');
addParameter(p,'excludenoise',true)
addParameter(p,'loadpc',false)
addParameter(p,'chmaskradius',55) % radius(um) within which the templates height are used to estimate position
addParameter(p,'getclufeature',true)
addParameter(p,'exportdir','')
parse(p,phydir,varargin{:});
phydir = p.Results.phydir;
excludenoise = p.Results.excludenoise;
loadpc = p.Results.loadpc;
chmaskradius = p.Results.chmaskradius;
getclufeature = p.Results.getclufeature;
exportdir = p.Results.exportdir;
%%
    function [cids, cgs] = readClusterGroupsCSV(filename)
        % cids is length nClusters, the cluster ID numbers
        % cgs is length nClusters, the "cluster group":
        % - 0 = noise
        % - 1 = good
        % - 2 = mua
        % - 3 = unsorted
        
        fid = fopen(filename);
        C = textscan(fid, '%s%s');
        fclose(fid);
        
        cids = cellfun(@str2num, C{1}(2:end), 'uni', false);
        ise = cellfun(@isempty, cids);
        cids = [cids{~ise}];
        
        isUns = cellfun(@(x)strcmp(x,'unsorted'),C{2}(2:end));
        isMUA = cellfun(@(x)strcmp(x,'mua'),C{2}(2:end));
        isGood = cellfun(@(x)strcmp(x,'good'),C{2}(2:end));
        
        cgs = zeros(size(cids));
        cgs(isGood(~ise)) = 1;
        cgs(isMUA(~ise)) = 2;
        cgs(isUns(~ise)) = 3;
        
        [~,isort] = sort(cids);
        cids = cids(isort);
        cgs = cgs(isort);
    end
%% load Phy data
dataset=[];
if ~exist(phydir,'dir')
    warning('No Phy Dir:    %s', phydir);
    return
end
disp(['Reading Phy Dir:    ',phydir,'    ...']);
spike = NeuroAnalysis.Phy.loadphyparam(fullfile(phydir, 'params.py'));
spike.fs = spike.sample_rate;
spike.datapath = spike.dat_path;
spike.nch = spike.n_channels_dat;
spike.imecindex = num2str(spike.imecindex);
spike = rmfield(spike,'sample_rate');
spike = rmfield(spike,'dat_path');
spike = rmfield(spike,'n_channels_dat');
if isfield(spike,'rawdat_path')
    spike.rawdatapath = spike.rawdat_path;
    spike.nchraw = spike.n_channels_rawdat;
    spike = rmfield(spike,'rawdat_path');
    spike = rmfield(spike,'n_channels_rawdat');
end

isbinfilerange = exist(fullfile(phydir, 'binfilerange.npy'),'file');
if isbinfilerange
    binfilerange = readNPY(fullfile(phydir, 'binfilerange.npy'));
end
if exist(fullfile(phydir, 'dshift.npy'),'file')
    spike.dshift = readNPY(fullfile(phydir, 'dshift.npy'));
    spike.yblk = readNPY(fullfile(phydir, 'yblk.npy'));
end

spikeTimes = readNPY(fullfile(phydir, 'spike_times.npy'));
spikeTemplates = readNPY(fullfile(phydir, 'spike_templates.npy'));

if exist(fullfile(phydir, 'spike_clusters.npy'),'file')
    clu = readNPY(fullfile(phydir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

tempScalingAmps = readNPY(fullfile(phydir, 'amplitudes.npy'));
temps = readNPY(fullfile(phydir, 'templates.npy'));

spike.chanmap = int64(readNPY(fullfile(phydir, 'channel_map.npy')))+1;
spike.channelposition = readNPY(fullfile(phydir, 'channel_positions.npy'));
if exist(fullfile(phydir, 'channel_map_raw.npy'),'file')
    spike.chanmapraw = int64(readNPY(fullfile(phydir, 'channel_map_raw.npy')))+1;
end

spike.whiteningmatrix = readNPY(fullfile(phydir, 'whitening_mat.npy'));
spike.whiteningmatrixinv = readNPY(fullfile(phydir, 'whitening_mat_inv.npy'));
if exist(fullfile(phydir, 'whitening_mat_raw.npy'),'file')
    spike.whiteningmatrixraw = readNPY(fullfile(phydir, 'whitening_mat_raw.npy'));
    spike.whiteningmatrixinvraw = readNPY(fullfile(phydir, 'whitening_mat_inv_raw.npy'));
end

if loadpc
    pcFeat = readNPY(fullfile(phydir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
    pcFeatInd = readNPY(fullfile(phydir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

cgsKSFile = '';
if exist(fullfile(phydir, 'cluster_KSLabel.csv'),'file')
    cgsKSFile = fullfile(phydir, 'cluster_KSLabel.csv');
end
if exist(fullfile(phydir, 'cluster_KSLabel.tsv'),'file')
    cgsKSFile = fullfile(phydir, 'cluster_KSLabel.tsv');
end
% every cluster labeled by Kilosort and Phy
iscgKS = ~isempty(cgsKSFile);
if iscgKS
    [cidsKS, cgsKS] = readClusterGroupsCSV(cgsKSFile); % cluster ids already sorted
else
    error('No KS cluster group labels.');
end

cgsFile = '';
if exist(fullfile(phydir, 'cluster_groups.csv'),'file')
    cgsFile = fullfile(phydir, 'cluster_groups.csv');
end
if exist(fullfile(phydir, 'cluster_group.tsv'),'file')
    cgsFile = fullfile(phydir, 'cluster_group.tsv');
end
% clusters labeled by Phy operator
iscg = ~isempty(cgsFile);
if iscg
    [cids, cgs] = readClusterGroupsCSV(cgsFile); % cluster ids already sorted
    
    % noise clusters are labeled by Phy operator
    if excludenoise
        noisecluidx = cgs==0;
        noisecluid = cids(noisecluidx);
        
        cgs = cgs(~noisecluidx);
        cids = cids(~noisecluidx);
        vsi = ~ismember(clu, noisecluid);
        if loadpc
            spike.pcFeat = pcFeat(vsi, :,:);
            spike.pcFeatInd = pcFeatInd(vsi,:);
        end
        spikeTimes = spikeTimes(vsi);
        spikeTemplates = spikeTemplates(vsi);
        tempScalingAmps = tempScalingAmps(vsi);
        clu = clu(vsi);
        
        % merge Phy labels into KS labels
        if iscgKS
            noisecluksidx = ismember(cidsKS, noisecluid);
            cgsKS = cgsKS(~noisecluksidx);
            cidsKS = cidsKS(~noisecluksidx);
            for i=1:length(cids)
                cgsKS(cids(i)==cidsKS) = cgs(i);
            end
        end
    end
end
disp(['Reading Phy Dir:    ',phydir,'    Done.']);
%% Prepare Spike
switch (spike.sort_from)
    case 'kilosort2'
        % unwhiten template to accurately get template position and spread
        [tempcoords,spikeAmps,tempAmps,templates_maxwaveform,templates_waveform_feature]...
            = NeuroAnalysis.Base.templatefeature(temps,spike.whiteningmatrixinv,spike.channelposition,chmaskradius,spikeTemplates,tempScalingAmps,spike.fs);
        spike.templates = temps; % nTemplates x nTimePoints x nChannels, used to do spike sorting
        spike.templatesposition = tempcoords; % position from template spatial spread
        spike.templatesamplitude = tempAmps; % mean amplitude of spikes from a unwhiten, scaled template
        spike.templateswaveform = templates_maxwaveform; % unwhiten waveform
        spike.templateswaveformfeature = templates_waveform_feature;
        
        spike.time = NeuroAnalysis.Base.sample2time(double(spikeTimes),spike.fs,spike.secondperunit);
        spike.template = int64(spikeTemplates)+1;
        spike.templatescale = tempScalingAmps;
        spike.amplitude = spikeAmps; % scaled template amplitude
        
        spike.cluster = int64(clu)+1;
        spike.clusterid = unique(spike.cluster,'sorted');
        spike.clustergood = cgsKS; % match the already sorted cluster ids
        
        % Cluster may not map 1:1 to template, so we extract cluster waveform from data and get waveform feature
        [~,fn,fe] = fileparts(spike.datapath);
        dirparts = strsplit(phydir,filesep);
        bindir = join(dirparts(1:end-1),filesep);
        binpath = fullfile(bindir{1},[fn,fe]);
        if exist(binpath,'file') && getclufeature
            mmapbinfile = memmapfile(binpath,'Format',{'int16',[spike.nch,spike.nsample],'ap'});
            % binfile is raw data, so use identity whiteningmatrixinv
            [cluwaveforms,clucoords,clumaxwaveform,cluwaveformfeature] = NeuroAnalysis.Base.clusterfeature(mmapbinfile,...
                double(spikeTimes),spike.cluster,spike.clusterid,spike.chanmap,eye(size(spike.whiteningmatrixinv)),spike.channelposition,chmaskradius,spike.fs);
            spike.clusterwaveforms = cluwaveforms; % mean waveform on channels
            spike.clusterposition = clucoords; % position from cluster spatial spread
            spike.clusterwaveform = clumaxwaveform; % mean waveform on max amplitude channel
            spike.clusterwaveformfeature = cluwaveformfeature;
        end
    case 'kilosort3'
        % use correct whiteningmatrixinv to restore unwhiten template waveform
        [tempcoords,spikeAmps,tempAmps,templates_maxwaveform,templates_waveform_feature]...
            = NeuroAnalysis.Base.templatefeature(temps,spike.whiteningmatrixinvraw,spike.channelposition,chmaskradius,spikeTemplates,tempScalingAmps,spike.fs);
        spike.templates = temps; % nTemplates x nTimePoints x nChannels, used to do spike sorting
        spike.templatesposition = tempcoords; % position from template spatial spread
        spike.templatesamplitude = tempAmps; % mean amplitude of spikes from a unwhiten, scaled template
        spike.templateswaveform = templates_maxwaveform; % unwhiten waveform
        spike.templateswaveformfeature = templates_waveform_feature;
        
        spike.time = NeuroAnalysis.Base.sample2time(double(spikeTimes),spike.fs,spike.secondperunit);
        spike.template = int64(spikeTemplates)+1;
        spike.templatescale = tempScalingAmps;
        spike.amplitude = spikeAmps; % scaled template amplitude
        
        spike.cluster = int64(clu)+1;
        spike.clusterid = unique(spike.cluster,'sorted');
        spike.clustergood = cgsKS; % match the already sorted cluster ids
        
        % Cluster may not map 1:1 to template, so we extract cluster waveform from data and get waveform feature
        if exist(spike.datapath,'file') && getclufeature
            mmapbinfile = memmapfile(spike.datapath,'Format',{'int16',[spike.nch,spike.nsample],'ap'});
            % binfile is whiten data, so use raw whiteningmatrixinv
            [cluwaveforms,clucoords,clumaxwaveform,cluwaveformfeature] = NeuroAnalysis.Base.clusterfeature(mmapbinfile,...
                double(spikeTimes),spike.cluster,spike.clusterid,spike.chanmap,spike.whiteningmatrixinvraw,spike.channelposition,chmaskradius,spike.fs);
            spike.clusterwaveforms = cluwaveforms; % mean waveform on channels
            spike.clusterposition = clucoords; % position from cluster spatial spread
            spike.clusterwaveform = clumaxwaveform; % mean waveform on max amplitude channel
            spike.clusterwaveformfeature = cluwaveformfeature;
        end
end
spike.qcversion='phy';
spike.qc = [];
%% Merge to Dataset
fieldtomerge = ['spike',spike.imecindex,'_',spike.sort_from];
dspath = split(spike.dataset_path,', ');
if length(dspath)>1
    % phy result is from concat binary file, need to split it for each dataset
    if isbinfilerange
        for i=1:length(dspath)
            merge2dataset(dspath{i},fieldtomerge,NeuroAnalysis.Base.splitspike(spike,binfilerange(i:i+1)));
        end
    else
        warning('Need `binfilerange` file to split result to each dataset.');
        return
    end
else
    merge2dataset(dspath{1},fieldtomerge,spike);
end
% data have been merged, no need to export anymore.
dataset.earlyfinish = true;
%%
    function merge2dataset(datasetpath,mergefield,spike)
        if ~exist(datasetpath,'file')
            warning('No Corresponding Dataset to Merge to:    %s', datasetpath);
            return
        end
        disp(['Merge Phy Result to Dataset:    ',datasetpath,'    ...']);
        odataset = matfile(datasetpath,'Writable',true);
        if isfield(odataset,mergefield)
            odataset.(mergefield) = NeuroAnalysis.Base.copyStructFields(spike,odataset.(mergefield));
        else
            odataset.(mergefield) = spike;
        end
        disp(['Merge Phy Result to Dataset:    ',datasetpath,'    Done.']);
    end
end