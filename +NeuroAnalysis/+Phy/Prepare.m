function [dataset] = Prepare(phykilosortdir,varargin)
%PREPARE Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'phykilosortdir');
addParameter(p,'ExcludeNoise',true)
addParameter(p,'LoadPC',false)
addParameter(p,'exportdir','')
parse(p,phykilosortdir,varargin{:});
phykilosortdir = p.Results.phykilosortdir;
excludenoise = p.Results.ExcludeNoise;
loadpc = p.Results.LoadPC;
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
        
        cgs(isGood) = 1;
        cgs(isMUA) = 2;
        cgs(isUns) = 3;
    end
%% load spike data
dataset=[];
if ~exist(phykilosortdir,'dir')
    warning('No Phy Dir:    %s', phykilosortdir);
    return
end
disp(['Reading Phy Dir:    ',phykilosortdir,'    ...']);
spike = NeuroAnalysis.Phy.loadphyparam(fullfile(phykilosortdir, 'params.py'));
spike.fs = spike.sample_rate;
spike = rmfield(spike,'sample_rate');

ss = readNPY(fullfile(phykilosortdir, 'spike_times.npy'));
spikeTemplates = readNPY(fullfile(phykilosortdir, 'spike_templates.npy'));

if exist(fullfile(phykilosortdir, 'spike_clusters.npy'),'file')
    clu = readNPY(fullfile(phykilosortdir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

tempScalingAmps = readNPY(fullfile(phykilosortdir, 'amplitudes.npy'));

if loadpc
    pcFeat = readNPY(fullfile(phykilosortdir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
    pcFeatInd = readNPY(fullfile(phykilosortdir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

cgsFile = '';
if exist(fullfile(phykilosortdir, 'cluster_groups.csv'),'file')
    cgsFile = fullfile(phykilosortdir, 'cluster_groups.csv');
end
if exist(fullfile(phykilosortdir, 'cluster_group.tsv'),'file')
    cgsFile = fullfile(phykilosortdir, 'cluster_group.tsv');
end
if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);
    
    if excludenoise
        noiseClusters = cids(cgs==0);
        
        ss = ss(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
        tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));
        
        if loadpc
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));
    end
else
    clu = spikeTemplates;
    
    cids = unique(spikeTemplates);
    cgs = 3*ones(size(cids));
end

coords = readNPY(fullfile(phykilosortdir, 'channel_positions.npy'));
temps = readNPY(fullfile(phykilosortdir, 'templates.npy'));
winv = readNPY(fullfile(phykilosortdir, 'whitening_mat_inv.npy'));
disp(['Reading Phy Dir:    ',phykilosortdir,'    Done.']);
%% Prepare Spike and Dataset
spike.time = double(ss)/spike.fs;
spike.template = int64(spikeTemplates)+1;
spike.cluster = int64(clu)+1;
spike.clusterid = int64(cids)+1;
spike.amplitude = tempScalingAmps;
spike.templates = temps;
spike.channelposition = coords;
spike.whiteningmatrixinv = winv;
spike.good = cgs;

spike.pcFeat = pcFeat;
spike.pcFeatInd = pcFeatInd;
spike.ischecked=true;

dataset.filepath = spike.dataset_path;
dataset.spike_kilosort = spike;
end