function [ dataset ] = KiloSort( dataset,ds )
%KILOSORT Summary of this function goes here
%   Detailed explanation goes here

if nargin==1
    ds='ap';
end

%% Concat binary files in time order into one binary file, then kilosort on it to get consistent cluster id across multiple binary files
if iscell(dataset)
    datasets = cellfun(@load,dataset,'uniformoutput',0);
    binfiles = cellfun(@(x)x.ap.meta.fileName,datasets,'uniformoutput',0);
    binfiledate = cellfun(@(x)x.ap.meta.fileDate,datasets);
    binfilensample = cellfun(@(x)x.ap.meta.nFileSamp,datasets);
    
    [~,isort]=sort(binfiledate);
    dataset = dataset(isort);
    datasets = datasets(isort);
    binfiles=binfiles(isort);
    binfilensample = binfilensample(isort);
    % make concat file name
    if strcmp(datasets{1}.ex.sourceformat, 'Stimulator')
        binrootdir = fileparts(binfiles{1});
        concatname=cell(1,length(binfiles)+1);
        concatname{1} = datasets{1}.source(1:8);  % Works for AE9, not AE4 yet!
        for i=1:length(binfiles)
            concatname{i+1} = datasets{i}.source(10:12);
        end
    else
        concatname=cell(size(binfiles));
        for i=1:length(binfiles)
            [binrootdir,fn,~] = fileparts(binfiles{i});
            tn = strsplit(fn,'_');
            tn(cellfun(@(x)isempty(x),tn))=[];
            concatname{i} = cellfun(@(x)x(1),tn(end-1:end));
        end
    end
    concatname=strjoin(concatname,'');
    concatfilepath = fullfile(binrootdir,concatname);
    if length(concatfilepath)>150
        concatfilepath = concatfilepath(1:150); % limit path name length on NTFS
    end
    concatfilepath = [concatfilepath,'.imec.ap.bin'];
    % demux channel groups
    nchsaved = datasets{1}.ap.meta.nSavedChans;
    dmxgroup = datasets{1}.ap.meta.dmxgroup;
    if datasets{1}.ap.meta.probeversion>0
        rmdc = false;
    else
        rmdc = true;
    end
    for g = 1:length(dmxgroup)
        dmxgroup{g} = setdiff(dmxgroup{g},datasets{1}.ap.meta.excludechans);
    end
    dmxgroup(cellfun(@(x)isempty(x),dmxgroup))=[];
    % concat binary files
    if exist(concatfilepath,'file')
        disp(['Use Existing Concat Binary File:    ',concatfilepath,'    ...']);
    else
        cfid=fopen(concatfilepath,'w');
        chunksample=1e7; % 1e7 Samples = 7.7GB for 385Chs, Int16
        for i=1:length(binfiles)
            fprintf('Concat Binary File:    %s    ...\n',binfiles{i});
            fid=fopen(binfiles{i},'r');
            remainingsample = binfilensample(i);
            while remainingsample>0
                if remainingsample - chunksample >=0
                    readsample = chunksample;
                else
                    readsample = remainingsample;
                end
                chunkdata = fread(fid,[nchsaved,readsample],'*int16');
                fprintf('        Demuxed CAR on chunk size:    %d    ...\n',readsample);
                chunkdata = NeuroAnalysis.Base.dmxcar(chunkdata,dmxgroup,rmdc);
                fwrite(cfid,chunkdata,'int16');
                remainingsample = remainingsample - readsample;
            end
            fclose(fid);
        end
        fclose(cfid);
        disp('Concat Binary Files        Done.');
    end
    % prepare a new dataset with concat binary file
    cdataset.secondperunit = datasets{1}.secondperunit;
    cdataset.filepath = dataset;
    cdataset.ap.meta = datasets{1}.ap.meta;
    cdataset.ap.meta.fileName=concatfilepath;
    cdataset.ap.meta.nFileSamp=sum(binfilensample);
    % demuxed CAR could help to remove very fast transient noise, and then normal CAR in kilosort could further reduce other noise.
    cdataset.car = 1;
    % time range [t(i), t(i+1)) for each file in the concat file
    cdataset.binfilerange = NeuroAnalysis.Base.sample2time(cumsum([1,binfilensample]),cdataset.ap.meta.fs,cdataset.secondperunit);
    clear datasets chunkdata % reclaim memory before KiloSort
    cdataset = NeuroAnalysis.SpikeGLX.KiloSort(cdataset);
    % split sorting result into each original dataset
    if isfield(cdataset,'spike_kilosort')
        disp('Split Sorting Result Into Dataset        ...');
        for i=1:length(dataset)
            odataset = matfile(dataset{i},'Writable',true);
            odataset.spike_kilosort = NeuroAnalysis.Base.splitspike(cdataset.spike_kilosort,cdataset.binfilerange(i:i+1));
        end
        disp('Split Sorting Result Into Dataset        Done.');
    end
    return;
end

%% Choose Kilosort Version
kilosortversion = '2';

rmpath(genpath('C:\Users\fff00\Kilosort-2.0'));
rmpath(genpath('C:\Users\fff00\Kilosort'));
switch kilosortversion
    case '2'
        addpath(genpath('C:\Users\fff00\Kilosort-2.0'));
    case '3'
        addpath(genpath('C:\Users\fff00\Kilosort'));
end

%% KiloSort ops
disp(['KiloSort ',kilosortversion,' Spike Sorting:    ',dataset.(ds).meta.fileName,'    ...']);

% the binary file
ops.fbinary = dataset.(ds).meta.fileName;

% the binary file folder
[binrootdir,binname,~] = fileparts(ops.fbinary);

% the probe channel map
[thisdir,~,~] = fileparts(mfilename('fullpath'));
ops.chanMap = neuropixelschmap(dataset.(ds).meta);

% exclude channels
excludechans = dataset.(ds).meta.excludechans;
if ~isempty(excludechans)
    disp(['Excluding Channels: ', num2str(excludechans)])
    ops.chanMap.connected(excludechans)=false;
end

% sample rate
ops.fs = dataset.(ds).meta.fs;

% frequency for high pass filtering
ops.fshigh = 300;

% frequency for low pass filtering
% narrow spike have significent power in high frequencies(up to 10kHz), if low pass
% filtering cuts off high freq component, the spike time would be shifted
% and spike shape be widen, so here we are not doing low pass.
% ops.fslow = 7000;

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 1/60;

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/60;

% number of samples to average over (annealed from first to second value)
ops.momentum = [20 400];

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30;

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8;

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [12 4];

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
ops.lam = 10;

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9;

% proc file on a fast SSD
ops.fproc = fullfile(binrootdir,['kilosort',kilosortversion,'_temp_wh.dat']);

% time range to sort
ops.trange = [0 Inf];

% total number of channels in your recording
ops.NchanTOT = double(dataset.(ds).meta.nSavedChans);

% common average referencing by median
ops.CAR = 1;
if isfield(dataset,'car')
    ops.CAR = dataset.car;
end

switch kilosortversion
    case '3'
        % spatial smoothness constant for registration
        ops.sig = 20;
        % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
        ops.nblocks = 5;
        
        ops.Th = [10 10];
        ops.lam = 20;
        ops.AUCsplit = 0.95;
end

%% danger, changing these settings can lead to fatal errors

% options for determining PCs
ops.spkTh               = -6;  % spike threshold in standard deviations (-6)
ops.reorder             = 1;   % whether to reorder batches for drift correction.
ops.nskip               = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64; % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 64*1024 + ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200; % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

%% run all the steps of kilosort
switch kilosortversion
    case '2'
        % preprocess data to create temp_wh.dat
        rez = preprocessDataSub(ops);
        % time-reordering as a function of drift
        rez = clusterSingleBatches(rez);
        % main tracking and template matching algorithm
        rez = learnAndSolve8b(rez);
        % OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
        % See issue 29: https://github.com/MouseLand/Kilosort2/issues/29
        % rez = remove_ks2_duplicate_spikes(rez,'channel_separation_um',60);
        % final merges
        rez = find_merges(rez, 1);
        % final splits by SVD
        rez = splitAllClusters(rez, 1);
        % final splits by amplitudes
        rez = splitAllClusters(rez, 0);
        % decide on cutoff
        rez = set_cutoff(rez);
    case '3'
        rez                = preprocessDataSub(ops);
        rez                = datashift2(rez, 1);
        [rez, st3, tF]     = extract_spikes(rez);
        rez                = template_learning(rez, tF, st3);
        [rez, st3, tF]     = trackAndSort(rez);
        rez                = final_clustering(rez, tF, st3);
        rez                = find_merges(rez, 1);
end

fprintf('Found %d / %d good units \n', sum(rez.good>0),length(rez.good));
disp(['KiloSort ',kilosortversion,' Spike Sorting:    ',dataset.(ds).meta.fileName,'    done.']);

%% Get sorted spikes
disp('Extract Spike Sorting Result    ...');
spike = [];%extractrez(rez,dataset.secondperunit,55);
disp('Extract Spike Sorting Result    done.');
if ~isempty(spike)
    dataset.spike_kilosort=spike;
end

phydir = fullfile(binrootdir,[binname,'_Kilosort',kilosortversion,'_Phy']);
if exist(phydir,'dir')
    rmdir(phydir,'s');
end
if ~exist(phydir,'dir')
    mkdir(phydir);
end
disp(['Save Spike Sorting Result for Phy in:    ',phydir,'    ...']);
switch kilosortversion
    case '2'
        rez2phy(rez,phydir,dataset,ds);
    case '3'
        rez2phy2(rez,phydir,dataset);
end
disp(['Save Spike Sorting Result for Phy in:    ',phydir,'    done.']);

%%
    function [chmap] = neuropixelschmap(meta)
        if meta.probeversion <= 1
            chmap.name = 'Neuropixels Phase3A/3B/1.0';
            chmap.connected = true(size(meta.roch));
            chmap.shankind = ones(size(chmap.connected));
            chmap.chanMap0ind = double(meta.robank*meta.acqApLfSy(1) + meta.roch);
            chmap.chanMap = chmap.chanMap0ind+1;
            dx = meta.probespacing(1);dy=meta.probespacing(2);
            cols = double(meta.savedcols);rows = double(meta.savedrows);
            
            if meta.probeversion <= 1 % Phase3A/3B/1.0 checkboard
                chmap.xcoords = arrayfun(@(r,c)(dx/2 + (c-1)*dx) - abs(mod(r,2)-1)*dx/2,rows,cols);
                chmap.ycoords = (rows-1)*dy;
            else % regular
                chmap.xcoords = (cols-1)*dx;
                chmap.ycoords = (rows-1)*dy;
            end
        end
    end
%% for kilosort 2
    function rez2phy(rez, savePath,dataset,ds)
        % pull out results from kilosort's rez to either return to workspace or to
        % save in the appropriate format for the phy GUI to run on. If you provide
        % a savePath it should be a folder, and you will need to have npy-matlab
        % available (https://github.com/kwikteam/npy-matlab)
        %
        % This is a modified version of the `rezToPhy` function in Kilosort
        
        if nargin==3
            ds='ap';
        end
        
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
        
        fs = dir(fullfile(savePath, '*.npy'));
        for i = 1:length(fs)
            delete(fullfile(savePath, fs(i).name));
        end
        if exist(fullfile(savePath, '.phy'), 'dir')
            rmdir(fullfile(savePath, '.phy'), 's');
        end
        
        % spikeTimes will be in samples, not seconds
        spikeTimes = uint64(rez.st3(:,1));
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
        U = rez.U;
        W = rez.W;
        nt0 = size(W,1);
        Nfilt = size(W,2);
        
        templates = zeros(Nchan, nt0, Nfilt, 'single');
        for iNN = 1:Nfilt
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
        
        if ~isempty(savePath)
            
            writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
            writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
            if size(rez.st3,2)>4
                writeNPY(uint32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
            else
                writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
            end
            writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy'));
            writeNPY(templates, fullfile(savePath, 'templates.npy'));
            writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
            
            writeNPY(int32(chanMap0ind), fullfile(savePath, 'channel_map.npy'));
            writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));
            
            writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
            writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
            writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
            writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
            
            writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
            writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
            
            if isfield(rez, 'simScore')
                similarTemplates = rez.simScore;
                writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
            end
            
            % save a list of "good" clusters for Phy
            fileID = fopen(fullfile(savePath, 'cluster_KSLabel.tsv'),'w');
            fprintf(fileID, 'cluster_id%sKSLabel', char(9));
            fprintf(fileID, char([13 10]));
            
            fileIDCP = fopen(fullfile(savePath, 'cluster_ContamPct.tsv'),'w');
            fprintf(fileIDCP, 'cluster_id%sContamPct', char(9));
            fprintf(fileIDCP, char([13 10]));
            
            fileIDA = fopen(fullfile(savePath, 'cluster_Amplitude.tsv'),'w');
            fprintf(fileIDA, 'cluster_id%sAmplitude', char(9));
            fprintf(fileIDA, char([13 10]));
            
            rez.est_contam_rate(isnan(rez.est_contam_rate)) = 1;
            for j = 1:length(rez.good)
                if rez.good(j)
                    fprintf(fileID, '%d%sgood', j-1, char(9));
                else
                    fprintf(fileID, '%d%smua', j-1, char(9));
                end
                fprintf(fileID, char([13 10]));
                
                fprintf(fileIDCP, '%d%s%.1f', j-1, char(9), rez.est_contam_rate(j)*100);
                fprintf(fileIDCP, char([13 10]));
                
                fprintf(fileIDA, '%d%s%.1f', j-1, char(9), tempAmps(j));
                fprintf(fileIDA, char([13 10]));
            end
            fclose(fileID);
            fclose(fileIDCP);
            fclose(fileIDA);
            
            % make params file
            if ~exist(fullfile(savePath,'params.py'),'file')
                fid = fopen(fullfile(savePath,'params.py'), 'w');
                
                [~, fname, ext] = fileparts(rez.ops.fbinary);
                
                % save dataset path so that phy sorting result could be loaded back into dataset
                if iscell(dataset.filepath)
                    datasetpath = strjoin(dataset.filepath,', ');
                else
                    datasetpath = dataset.filepath;
                end
                fprintf(fid,['dataset_path = ''',replace(datasetpath,'\','\\'),'''\n']);
                fprintf(fid,'secondperunit = %.32f\n',dataset.secondperunit);
                if isfield(dataset,'binfilerange')
                    writeNPY(dataset.binfilerange, fullfile(savePath, 'binfilerange.npy'));
                end
                fprintf(fid,['dat_path = ''','../',fname ext '''\n']); % phy folder in the sam folder of binaries
                fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
                fprintf(fid,'nsample = %i\n',dataset.(ds).meta.nFileSamp);
                fprintf(fid,'dtype = ''int16''\n');
                fprintf(fid,'offset = 0\n');
                fprintf(fid,'sample_rate = %.32f\n',rez.ops.fs);
                fprintf(fid,['sort_from = ''Kilosort_v',kilosortversion,'''\n']);
                if isfield(rez.ops,'fshigh')
                    hp='True';
                else
                    hp='False';
                end
                fprintf(fid,['hp_filtered = ', hp]);
                fclose(fid);
            end
        end
    end
%% for kilosort 3
    function rez2phy2(rez, savePath,dataset)
        % pull out results from kilosort's rez to either return to workspace or to
        % save in the appropriate format for the phy GUI to run on. If you provide
        % a savePath it should be a folder, and you will need to have npy-matlab
        % available (https://github.com/kwikteam/npy-matlab)
        %
        % This is a modified version of the `rezToPhy2` function in Kilosort
        
        [~, Nfilt, Nrank] = size(rez.W);
        % for Phy, we need to pad the spikes with zeros so the spikes are aligned to the center of the window
        rez.Wphy = cat(1, zeros(1+rez.ops.nt0min, Nfilt, Nrank), rez.W);
        
        rez.W = gather(single(rez.Wphy));
        rez.U = gather(single(rez.U));
        rez.mu = gather(single(rez.mu));
        
        if size(rez.st3,2)>4
            rez.st3 = rez.st3(:,1:4);
        end
        
        [~, isort]   = sort(rez.st3(:,1), 'ascend');
        rez.st3      = rez.st3(isort, :);
        if ~isempty(rez.cProj)
            rez.cProj    = rez.cProj(isort, :);
            rez.cProjPC  = rez.cProjPC(isort, :, :);
        end
        
        fs = dir(fullfile(savePath, '*.npy'));
        for i = 1:length(fs)
            delete(fullfile(savePath, fs(i).name));
        end
        if exist(fullfile(savePath, '.phy'), 'dir')
            rmdir(fullfile(savePath, '.phy'), 's');
        end
        
        % spikeTimes will be in samples, not seconds
        spikeTimes = uint64(rez.st3(:,1));
        % account for ops.trange(1) to accomodate real time
        spikeTimes = spikeTimes - rez.ops.trange(1)*rez.ops.fs;
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
        U = rez.U;
        W = rez.W;
        nt0 = size(W,1);
        Nfilt = size(W,2);
        
        templates = zeros(Nchan, nt0, Nfilt, 'single');
        for iNN = 1:Nfilt
            templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
        end
        templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
        templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial
        
        templateFeatures = rez.cProj;
        templateFeatureInds = uint32(rez.iNeigh);
        pcFeatures = rez.cProjPC;
        pcFeatureInds = uint32(rez.iNeighPC);
        
        %         whiteningMatrix = rez.Wrot/rez.ops.scaleproc;
        whiteningMatrix = eye(size(rez.Wrot)) / rez.ops.scaleproc;
        whiteningMatrixInv = whiteningMatrix^-1;
        
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
        tempAmps = zeros(numel(rez.mu),1);
        tempAmps(tids) = ta; % because ta only has entries for templates that had at least one spike
        gain = getOr(rez.ops, 'gain', 1);
        tempAmps = gain*tempAmps'; % for consistency, make first dimension template number
        
        if ~isempty(savePath)
            
            writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
            writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
            if size(rez.st3,2)>4
                writeNPY(uint32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
            else
                writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
            end
            writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy'));
            writeNPY(templates, fullfile(savePath, 'templates.npy'));
            writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
            
            writeNPY(int32(chanMap0ind), fullfile(savePath, 'channel_map.npy'));
            writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));
            
            if ~isempty(templateFeatures)
                writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
                writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
                writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
                writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
            end
            
            writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
            writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
            
            if isfield(rez, 'simScore')
                similarTemplates = rez.simScore;
                writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
            end
            
            % save a list of "good" clusters for Phy
            fileID = fopen(fullfile(savePath, 'cluster_KSLabel.tsv'),'w');
            fprintf(fileID, 'cluster_id%sKSLabel', char(9));
            fprintf(fileID, char([13 10]));
            
            fileIDCP = fopen(fullfile(savePath, 'cluster_ContamPct.tsv'),'w');
            fprintf(fileIDCP, 'cluster_id%sContamPct', char(9));
            fprintf(fileIDCP, char([13 10]));
            
            fileIDA = fopen(fullfile(savePath, 'cluster_Amplitude.tsv'),'w');
            fprintf(fileIDA, 'cluster_id%sAmplitude', char(9));
            fprintf(fileIDA, char([13 10]));
            
            rez.est_contam_rate(isnan(rez.est_contam_rate)) = 1;
            for j = 1:length(rez.good)
                if rez.good(j)
                    fprintf(fileID, '%d%sgood', j-1, char(9));
                else
                    fprintf(fileID, '%d%smua', j-1, char(9));
                end
                fprintf(fileID, char([13 10]));
                
                fprintf(fileIDCP, '%d%s%.1f', j-1, char(9), rez.est_contam_rate(j)*100);
                fprintf(fileIDCP, char([13 10]));
                
                fprintf(fileIDA, '%d%s%.1f', j-1, char(9), tempAmps(j));
                fprintf(fileIDA, char([13 10]));
            end
            fclose(fileID);
            fclose(fileIDCP);
            fclose(fileIDA);
            
            % Duplicate "KSLabel" as "group", a special metadata ID for Phy, so that
            % filtering works as expected in the cluster view
            KSLabelFilename = fullfile(savePath, 'cluster_KSLabel.tsv');
            copyfile(KSLabelFilename, fullfile(savePath, 'cluster_group.tsv'));
            
            % make params file
            if ~exist(fullfile(savePath,'params.py'),'file')
                fid = fopen(fullfile(savePath,'params.py'), 'w');
                
                [~, fname, ext] = fileparts(rez.ops.fbinary);
                
                % save dataset path so that phy sorting result could be loaded back into dataset
                if iscell(dataset.filepath)
                    datasetpath = strjoin(dataset.filepath,', ');
                else
                    datasetpath = dataset.filepath;
                end
                fprintf(fid,['dataset_path = ''',replace(datasetpath,'\','\\'),'''\n']);
                fprintf(fid,'secondperunit = %.32f\n',dataset.secondperunit);
                if isfield(dataset,'binfilerange')
                    writeNPY(dataset.binfilerange, fullfile(savePath, 'binfilerange.npy'));
                end
                fprintf(fid,['dat_path = ''', strrep(rez.ops.fproc, '\', '/') '''\n']);
                %                 fprintf(fid,['dat_path = ''','../',fname ext '''\n']); % phy folder usually in the sam folder of binaries
                fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
                fprintf(fid,'nsample = %i\n',dataset.ap.meta.nFileSamp);
                fprintf(fid,'dtype = ''int16''\n');
                fprintf(fid,'offset = 0\n');
                fprintf(fid,'sample_rate = %.32f\n',rez.ops.fs);
                fprintf(fid,['sort_from = ''Kilosort_v',kilosortversion,'''\n']);
                if isfield(rez.ops,'fshigh')
                    hp='True';
                else
                    hp='False';
                end
                fprintf(fid,['hp_filtered = ', hp]);
                fclose(fid);
            end
        end
    end
%% extract spike directly from kilosort result
% chmaskradius(um) within which the templates height are used to estimate position
    function [spike]=extractrez(rez,secondperunit,chmaskradius)
        W = gather(single(rez.Wphy));
        U = gather(single(rez.U));
        nt0 = size(W,1);
        Nfilt = size(W,2);
        Nch = rez.ops.Nchan;
        
        % each template in spatial(electrodes) and temporal(waveform) dimention
        temps = zeros(Nch, nt0, Nfilt, 'single');
        for iNN = 1:Nfilt
            temps(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
        end
        temps = permute(temps, [3 2 1]); % nTemplates x nSamples x nChannels
        
        w = rez.Wrot/rez.ops.scaleproc;
        winv = w^-1;
        coords = [rez.xcoords(:) rez.ycoords(:)];
        spikeTimes = rez.st3(:,1);
        spikeTemplates = rez.st3(:,2);
        tempScalingAmps = rez.st3(:,3);
        chanMap = rez.ops.chanMap(:);
        
        rez.ops.igood = gather(rez.ops.igood);
        spike.ops = rez.ops;
        spike.fs = rez.ops.fs;
        
        [tempcoords,spikeAmps,tempAmps,templates_maxwaveform,templates_waveform_feature]...
            = NeuroAnalysis.Base.templatefeature(temps,winv,coords,chmaskradius,spikeTemplates-1,tempScalingAmps,spike.fs);
        % Templates feature
        spike.templates = temps;
        spike.templatesposition = tempcoords;
        spike.templatesamplitude = tempAmps;
        spike.templateswaveform = templates_maxwaveform;
        spike.templateswaveformfeature = templates_waveform_feature;
        
        spike.chanmap = int32(chanMap);
        spike.channelposition = coords;
        spike.whiteningmatrix = w;
        spike.whiteningmatrixinv = winv;
        
        % times for each spike, t0 is the first sample in the data stream
        spike.time = NeuroAnalysis.Base.sample2time(spikeTimes,spike.fs,secondperunit);
        % ids of template on which each spike is extracted
        spike.template = int64(spikeTemplates);
        % template scaling for each spike
        spike.templatescale = tempScalingAmps;
        % scaled unwhiten template amplitude for each spike
        spike.amplitude = spikeAmps;
        
        % ids of cluster where each is expected to be a single cell, here
        % regarding each unique template as a single cell
        spike.cluster = spike.template;
        % unique cluser ids(already sorted in `unique`)
        spike.clusterid = unique(spike.cluster);
        spike.clustergood = rez.good;
        
        spike.qcversion = ['Kilosort_v',kilosortversion];
        spike.qc = [];
    end
end