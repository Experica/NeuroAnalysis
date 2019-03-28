function [ dataset ] = KiloSort( dataset )
%KILOSORT Summary of this function goes here
%   Detailed explanation goes here

%% KiloSort2 ops
[thisdir,~,~] = fileparts(mfilename('fullpath'));
ops.chanMap             = fullfile(thisdir,'neuropixPhase3B2_kilosortChanMap.mat');
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if no chanMap file

% sample rate
ops.fs = 30000;

% frequency for high pass filtering (150)
ops.fshigh = 150;

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.1;

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4];

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
ops.lam = 10;

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9;

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/50;

% number of samples to average over (annealed from first to second value)
ops.momentum = [20 400];

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30;

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8;

% the binary file
ops.fbinary = dataset.spike.meta.fileName;
% the binary file folder
[rootdir,filename,ext] = fileparts(ops.fbinary);

ops.fproc       = 'E:\temp_wh.dat'; % proc file on a fast SSD
ops.trange = [0 Inf]; % time range to sort
ops.NchanTOT    = 385; % total number of channels in your recording
ops.CAR = 1;
%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction.
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available
%% this block runs all the steps of the algorithm

% is there a channel map file in this folder?
fs = dir(fullfile(rootdir, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootdir, fs(1).name);
end

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);
%save(fullfile(rootdir, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

%% Get Sorted spikes
% write to Phy
fprintf('Saving results to Phy  \n')
%rezToPhy(rez, rootdir);

spike = extractrez(rez);
if ~isempty(spike)
    dataset.spike=spike;
end

    function [spike]=extractrez(rez)
        % spikeTimes will be in samples, not seconds
        rez.W = gather(single(rez.Wphy));
        rez.U = gather(single(rez.U));
        rez.mu = gather(single(rez.mu));
        
        if size(rez.st3,2)>4
            rez.st3 = rez.st3(:,1:4);
        end
        
        [~, isort]   = sort(rez.st3(:,1), 'ascend');
        rez.st3      = rez.st3(isort, :);
        rez.cProj    = rez.cProj(isort, :);
        rez.cProjPC  = rez.cProjPC(isort, :, :);
        
        spikeTimes = uint64(rez.st3(:,1));
        % [spikeTimes, ii] = sort(spikeTimes);
        spikeTemplates = uint32(rez.st3(:,2));
        if size(rez.st3,2)>4
            spikeClusters = uint32(1+rez.st3(:,5));
        end
        amplitudes = rez.st3(:,3);
        
        Nchan = rez.ops.Nchan;
        
        xcoords     = rez.xcoords(:);
        ycoords     = rez.ycoords(:);
        chanMap     = rez.ops.chanMap(:);
        chanMap0ind = chanMap - 1;
        
        nt0 = size(rez.W,1);
        U = rez.U;
        W = rez.W;
        
        Nfilt = size(W,2);
        
        templates = zeros(Nchan, nt0, Nfilt, 'single');
        for iNN = 1:size(templates,3)
            templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
        end
        templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
        templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial
        
        templateFeatures = rez.cProj;
        templateFeatureInds = uint32(rez.iNeigh);
        pcFeatures = rez.cProjPC;
        pcFeatureInds = uint32(rez.iNeighPC);
        
        whiteningMatrix = rez.Wrot/rez.ops.scaleproc;
        whiteningMatrixInv = whiteningMatrix^-1;
        
        % here we compute the amplitude of every template...
        
        % unwhiten all the templates
        tempsUnW = zeros(size(templates));
        for t = 1:size(templates,1)
            tempsUnW(t,:,:) = squeeze(templates(t,:,:))*whiteningMatrixInv;
        end
        
        % The amplitude on each channel is the positive peak minus the negative
        tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));
        
        % The template amplitude is the amplitude of its largest channel
        tempAmpsUnscaled = max(tempChanAmps,[],2);
        
        % assign all spikes the amplitude of their template multiplied by their
        % scaling amplitudes
        spikeAmps = tempAmpsUnscaled(spikeTemplates).*amplitudes;
        
        % take the average of all spike amps to get actual template amps (since
        % tempScalingAmps are equal mean for all templates)
        ta = clusterAverage(spikeTemplates, spikeAmps);
        tids = unique(spikeTemplates);
        tempAmps(tids) = ta; % because ta only has entries for templates that had at least one spike
        gain = getOr(rez.ops, 'gain', 1);
        tempAmps = gain*tempAmps'; % for consistency, make first dimension template number
        
        
        spike.time = spikeTimes;
        spike.templates = uint32(spikeTemplates);
        
        
        if size(rez.st3,2)>4
            spike.clusters = uint32(spikeClusters);
        else
            spike.clusters = spike.templates;
        end
        spike.amplitudes = amplitudes;
        spike.templates_ind = templatesInds;
        
        chanMap0ind = int32(chanMap0ind);
        
        spike.channel_map = chanMap0ind;
        spike.channel_positions = [xcoords ycoords];
        
        spike.template_features = templateFeatures;
        spike.template_feature_ind = templateFeatureInds';
        spike.pc_features = pcFeatures;
        spike.pc_feature_ind = pcFeatureInds';
        
        spike.whitening_mat = whiteningMatrix;
        spike.whitening_mat_inv = whiteningMatrixInv;
        
        if isfield(rez, 'simScore')
            spike.similar_templates = rez.simScore;
        end
        spike.sortmethod = 'KiloSort2';
        spike.sortparams = rez.ops;
        
    end
end