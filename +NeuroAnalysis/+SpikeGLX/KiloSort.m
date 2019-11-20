function [ dataset ] = KiloSort( dataset )
%KILOSORT Summary of this function goes here
%   Detailed explanation goes here

%% Concat binary files in time order to get consistent cluster id
if iscell(dataset)
    datasets = cellfun(@load,dataset,'uniformoutput',0);
    binfiles = cellfun(@(x)x.ap.meta.fileName,datasets,'uniformoutput',0);
    binfiledate = cellfun(@(x)x.ap.meta.fileDate,datasets);
    binfilensample= cellfun(@(x)x.ap.meta.nFileSamp,datasets);
    
    [~,isort]=sort(binfiledate);
    binfiles=binfiles(isort);
    binfilensample = binfilensample(isort);
    dataset = dataset(isort);
    datasets = datasets(isort);
    % get concat file name
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
            [binrootdir,concatname{i},~] = fileparts(binfiles{i});
        end
    end
    concatname=strjoin(concatname,'__');
    concatfilepath = fullfile(binrootdir,concatname);
    if length(concatfilepath)>240
        concatfilepath = concatfilepath(1:240); % limit path name length on NTFS
    end
    concatfilepath = [concatfilepath,'.bin'];
    % concat binary files
    nch = datasets{1}.ap.meta.nSavedChans;
    if exist(concatfilepath,'file')
        disp(['Use Existing Concat Binary File:    ',concatfilepath,'    ...']);
    else
        cfid=fopen(concatfilepath,'w');
        chucksample=1000000; % 1MSample ~ 1GB for 500Chs
        for i=1:length(binfiles)
            fprintf('Concat Binary File:    %s    ...\n',binfiles{i});
            fid=fopen(binfiles{i},'r');
            remainingsample = binfilensample(i);
            while remainingsample>0
                if remainingsample - chucksample >=0
                    readsample = chucksample;
                else
                    readsample = remainingsample;
                end
                readdata = fread(fid,[nch,readsample],'*int16');
                fwrite(cfid,readdata,'int16');
                remainingsample = remainingsample - readsample;
            end
            fclose(fid);
        end
        fclose(cfid);
        disp('Concat Binary Files        Done.');
    end
    % concat binary file kilosort
    cdataset.secondperunit = datasets{1}.secondperunit;
    cdataset.filepath = dataset;
    cdataset.ap.meta = datasets{1}.ap.meta;
    cdataset.ap.meta.fileName=concatfilepath;
    cdataset.binfilerange = NeuroAnalysis.Base.sample2time(cumsum([1,binfilensample]),cdataset.ap.meta.fs,cdataset.secondperunit);
    clear datasets % reclaim memory before KiloSort
    cdataset = NeuroAnalysis.SpikeGLX.KiloSort(cdataset);
    % split sorting result into each original dataset
    disp('Split Sorting Result Into Dataset        ...');
    spike = cdataset.spike_kilosort;
    for i=1:length(dataset)
        odataset = matfile(dataset{i},'Writable',true);
        odataset.spike_kilosort = NeuroAnalysis.Base.splitspike(spike,cdataset.binfilerange(i:i+1));
    end
    disp('Split Sorting Result Into Dataset        Done.');
    return;
end

%% KiloSort2 ops
disp(['KiloSort2 Spike Sorting:    ',dataset.ap.meta.fileName,'    ...']);

% the binary file
ops.fbinary = dataset.ap.meta.fileName;

% the binary file folder
[binrootdir,binname,~] = fileparts(ops.fbinary);

% the probe channel map
[thisdir,~,~] = fileparts(mfilename('fullpath'));
ops.chanMap = load(fullfile(thisdir,'NeuropixelPhase3A_KilosortChannelMap.mat'));

% bad channels
if dataset.ap.meta.rorefch(1)==0
    badch = dataset.ap.meta.refch;
else
    badch=dataset.ap.meta.rorefch(1);
end
if isfield(dataset.ap.meta,'badch')
    badch = union(badch,dataset.ap.meta.badch);
end
if ~isempty(badch)
    disp(['Excluding Bad Channels: ', num2str(badch)])
    ops.chanMap.connected(badch)=0;
end

% sample rate
ops.fs = dataset.ap.meta.fs;

% frequency for high pass filtering
ops.fshigh = 300;

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.01;

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4];

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
ops.lam = 10;

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.95;

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/100;

% number of samples to average over (annealed from first to second value)
ops.momentum = [20 400];

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30;

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8;

% proc file on a fast SSD
ops.fproc = fullfile(binrootdir,'kilosort_temp_wh.dat');

% time range to sort
ops.trange = [0 Inf];

% total number of channels in your recording
ops.NchanTOT = double(dataset.ap.meta.nSavedChans);

% common average referencing by median
ops.CAR = 1;

%% danger, changing these settings can lead to fatal errors

% options for determining PCs
ops.spkTh           = -6;  % spike threshold in standard deviations (-6)
ops.reorder         = 1;   % whether to reorder batches for drift correction.
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64; % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 64*1024 + ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200; % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available
%% this block runs all the steps of the algorithm

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

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

fprintf('found %d / %d good units \n', sum(rez.good>0),length(rez.good))

%% Get sorted spikes
spike = extractrez(rez,dataset.secondperunit);
if ~isempty(spike)
    dataset.spike_kilosort=spike;
end

phydir = fullfile(binrootdir,[binname,'_Phy']);
if exist(phydir,'dir')
    rmdir(phydir,'s');
end
if ~exist(phydir,'dir')
    mkdir(phydir);
end
rez2phy(rez,phydir,dataset);

if exist(ops.fproc, 'file')
    delete(ops.fproc);
end
disp(['KiloSort2 Spike Sorting:    ',dataset.ap.meta.fileName,'    done.']);
%%
    function rez2phy(rez, savePath,dataset)
        % pull out results from kilosort's rez to either return to workspace or to
        % save in the appropriate format for the phy GUI to run on. If you provide
        % a savePath it should be a folder, and you will need to have npy-matlab
        % available (https://github.com/kwikteam/npy-matlab)
        
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
        
        % ix = rez.st3(:,4)>12;
        % rez.st3 = rez.st3(ix, :);
        % rez.cProj = rez.cProj(ix, :);
        % rez.cProjPC = rez.cProjPC(ix, :,:);
        
        fs = dir(fullfile(savePath, '*.npy'));
        for i = 1:length(fs)
            delete(fullfile(savePath, fs(i).name));
        end
        if exist(fullfile(savePath, '.phy'), 'dir')
            rmdir(fullfile(savePath, '.phy'), 's');
        end
        
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
            
            chanMap0ind = int32(chanMap0ind);
            
            writeNPY(chanMap0ind, fullfile(savePath, 'channel_map.npy'));
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
                fprintf(fid,['dat_path = ''','../',fname ext '''\n']); % phy folder usually in the sam folder of binaries
                fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
                fprintf(fid,'dtype = ''int16''\n');
                fprintf(fid,'offset = 0\n');
                if mod(rez.ops.fs,1)
                    fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
                else
                    fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
                end
                fprintf(fid,'hp_filtered = False');
                fclose(fid);
            end
        end
    end
%%
    function [spike]=extractrez(rez,secondperunit)
        spike.fs = rez.ops.fs;
        spike.time = NeuroAnalysis.Base.sample2time(rez.st3(:,1),spike.fs,secondperunit);
        spike.template = int64(rez.st3(:,2));
        spike.cluster = spike.template;
        spike.clusterid = unique(spike.cluster); % already sorted, match nTemplates
        spike.amplitude = rez.st3(:,3);
        
        rez.W = gather(single(rez.Wphy));
        rez.U = gather(single(rez.U));
        nt0 = size(rez.W,1);
        U = rez.U;
        W = rez.W;
        
        Nch = rez.ops.Nchan;
        Nfilt = size(W,2);
        templates = zeros(Nch, nt0, Nfilt, 'single');
        for iNN = 1:size(templates,3)
            templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
        end
        spike.clustertemplates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
        spike.chanmap = int64(rez.ops.chanMap(:));
        spike.channelposition = [rez.xcoords(:) rez.ycoords(:)];
        spike.whiteningmatrix = rez.Wrot/rez.ops.scaleproc;
        spike.whiteningmatrixinv = spike.whiteningmatrix^-1;
        spike.clustergood = rez.good;
        
        rez.ops.igood = gather(rez.ops.igood);
        spike.ops = rez.ops;
        spike.ischecked = false;
    end
end