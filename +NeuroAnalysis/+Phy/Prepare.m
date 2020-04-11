function [dataset] = Prepare(phydir,varargin)
%PREPARE Read Phy result in phydir and merge back to dataset
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'phydir');
addParameter(p,'excludenoise',true)
addParameter(p,'loadpc',false)
addParameter(p,'chmaskradius',65) % radius(um) within which the templates height are used to estimate position
addParameter(p,'getclufeature',false)
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
        
        cgs(isGood) = 1;
        cgs(isMUA) = 2;
        cgs(isUns) = 3;
        
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
spike.nch = spike.n_channels_dat;
spike = rmfield(spike,'sample_rate');
spike = rmfield(spike,'n_channels_dat');

isbinfilerange = exist(fullfile(phydir, 'binfilerange.npy'),'file');
if isbinfilerange
    binfilerange = readNPY(fullfile(phydir, 'binfilerange.npy'));
end

ss = readNPY(fullfile(phydir, 'spike_times.npy'));
spikeTemplates = readNPY(fullfile(phydir, 'spike_templates.npy'));

if exist(fullfile(phydir, 'spike_clusters.npy'),'file')
    clu = readNPY(fullfile(phydir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

tempScalingAmps = readNPY(fullfile(phydir, 'amplitudes.npy'));
chmap = readNPY(fullfile(phydir, 'channel_map.npy'));
coords = readNPY(fullfile(phydir, 'channel_positions.npy'));
temps = readNPY(fullfile(phydir, 'templates.npy'));
w = readNPY(fullfile(phydir, 'whitening_mat.npy'));
winv = readNPY(fullfile(phydir, 'whitening_mat_inv.npy'));

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
        vi = ~ismember(clu, noisecluid);
        if loadpc
            pcFeat = pcFeat(vi, :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end
        ss = ss(vi);
        spikeTemplates = spikeTemplates(vi);
        tempScalingAmps = tempScalingAmps(vi);
        clu = clu(vi);
        
        % overwrite KS labels with Phy labels
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
switch (spike.from)
    case 'kilosort'
        [tempcoords,spikeAmps,tempAmps,templates_maxwaveform_chidx,templates_maxwaveform,templates_maxwaveform_feature]...
            = NeuroAnalysis.Base.templatefeature(temps,winv,coords,chmaskradius,spikeTemplates,tempScalingAmps,spike.fs);
        % Templates feature
        spike.templates = temps; % nTemplates x nTimePoints x nChannels
        spike.templatesposition = tempcoords;
        spike.templatesamplitude = tempAmps;
        spike.templateswaveform = templates_maxwaveform;
        spike.templateswaveformfeature = templates_maxwaveform_feature;
        
        spike.chanmap = chmap+1;
        spike.channelposition = coords;
        spike.whiteningmatrix = w;
        spike.whiteningmatrixinv = winv;
        spike.pcFeat = pcFeat;
        spike.pcFeatInd = pcFeatInd;
        
        spike.time = NeuroAnalysis.Base.sample2time(double(ss),spike.fs,spike.secondperunit);
        spike.template = int64(spikeTemplates)+1;
        spike.templatescale = tempScalingAmps;
        spike.amplitude = spikeAmps;
        
        spike.cluster = int64(clu)+1;
        spike.clusterid = unique(spike.cluster); % already sorted in `unique`
        spike.clustergood = cgsKS; % match the already sorted cluster ids
        
        % Cluster feature
        [~,fn,fe] = fileparts(spike.dat_path);
        dirparts = strsplit(phydir,filesep);
        bindir = join(dirparts(1:end-1),filesep);
        binpath = fullfile(bindir{1},[fn,fe]);
        if exist(binpath,'file') && getclufeature
            mpbinfile = memmapfile(binpath,'Format',{'int16',[spike.nch,spike.nsample],'ap'});
            [cluwaveform,cluwaveformfeature] = NeuroAnalysis.Base.clusterfeature(mpbinfile,...
                double(ss),spike.cluster,spike.clusterid,spike.template,chmap(templates_maxwaveform_chidx),spike.fs);
            spike.clusterwaveform = cluwaveform;
            spike.clusterwaveformfeature = cluwaveformfeature;
        end
end
spike.qcversion='Phy';
spike.qc = [];
%% Merge to Dataset
fieldtomerge = ['spike_',spike.from];
ds = split(spike.dataset_path,', ');
if length(ds)>1
    % result is from concat binary file, need to split it for each dataset
    if isbinfilerange
        for i=1:length(ds)
            merge2dataset(ds{i},fieldtomerge,NeuroAnalysis.Base.splitspike(spike,binfilerange(i:i+1)));
        end
    else
        warning('Need `binfilerange` file to split result to each dataset.');
        return
    end
else
    merge2dataset(ds{1},fieldtomerge,spike);
end
dataset.earlyfinish = true;
%%
    function merge2dataset(datasetpath,mergefield,spike)
        if ~exist(datasetpath,'file')
            warning('No Corresponding Dataset to Merge to:    %s', datasetpath);
            return
        end
        disp(['Merge Phy Result to Dataset:    ',datasetpath,'    ...']);
        odataset = matfile(datasetpath,'Writable',true);
        odataset.(mergefield) = NeuroAnalysis.Base.copyStructFields(spike,odataset.(mergefield));
        disp(['Merge Phy Result to Dataset:    ',datasetpath,'    Done.']);
    end
end