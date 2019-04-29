function [ dataset ] = KiloSort( dataset )
%KILOSORT Summary of this function goes here
%   Detailed explanation goes here

%% Concat binary files in time order to get consistent sorting
if iscell(dataset)
    datasets = cellfun(@load,dataset,'uniformoutput',0);
    binfiles = cellfun(@(x)x.ap.meta.fileName,datasets,'uniformoutput',0);
    binfiledate = cellfun(@(x)x.ap.meta.fileDate,datasets);
    binfilensample= cellfun(@(x)x.ap.meta.nFileSamp,datasets);
    nch = datasets{1}.ap.meta.nSavedChans;
    [~,isort]=sort(binfiledate);
    binfiles=binfiles(isort);
    binfilensample = binfilensample(isort);
    datasets=datasets(isort);
    % get concat file name
    concatname=cell(size(binfiles));
    for i=1:length(binfiles)
        [binrootdir,concatname{i},~] = fileparts(binfiles{i});
    end
    concatname=[strjoin(concatname,'__'),'.bin'];
    if length(concatname)>240
        concatname = concatname(1:240); % max file/folder name length on NTFS
    end
    concatfilepath = fullfile(binrootdir,concatname);
    % concat binary files
    if exist(concatfilepath,'file')
        disp(['Use Existing Concat Binary File:    ',concatfilepath,'    ...']);
    else
        cfid=fopen(concatfilepath,'w');
        chucksize=1000000; % 1MSample ~ 1GB for 500Chs
        for i=1:length(binfiles)
            fprintf('Concat Binary File:    %s    ...\n',binfiles{i});
            fid=fopen(binfiles{i},'r');
            for j=1:ceil(binfilensample(i)/chucksize)
                chuckdata=fread(fid,[nch,chucksize],'*int16');
                fwrite(cfid,chuckdata,'int16');
            end
            fclose(fid);
        end
        fclose(cfid);
        disp('Concat Binary Files        Done.');
    end
    % concat binary file kilosort
    dataset=struct;
    dataset.secondperunit = datasets{1}.secondperunit;
    dataset.ap.meta.fileName=concatfilepath;
    dataset.ap.meta.imSampRate = datasets{1}.ap.meta.imSampRate;
    dataset.ap.meta.nSavedChans= nch;
    dataset = NeuroAnalysis.SpikeGLX.KiloSort(dataset);
    % split sorting result into each original dataset
    disp('Split Sorting Result Into Dataset Files        ...');
    spike = dataset.spike_kilosort;
    binfilerange=NeuroAnalysis.Base.sample2time(cumsum([1,binfilensample]),spike.fs,dataset.secondperunit);
    for i=1:length(datasets)
        odataset = datasets{i};
        odataset.spike_kilosort = splitspike(spike,binfilerange(i:i+1));
        save(odataset.filepath,'-struct','odataset','-v7.3');
    end
    disp('Split Sorting Result Into Dataset Files        Done.');
    return;
end

disp(['KiloSort2 Spike Sorting:    ',dataset.ap.meta.fileName,'    ...']);
%% KiloSort2 ops
% the binary file
ops.fbinary = dataset.ap.meta.fileName;

% the binary file folder
[binrootdir,binname,~] = fileparts(ops.fbinary);

% the probe channel map
[thisdir,~,~] = fileparts(mfilename('fullpath'));
ops.chanMap = fullfile(thisdir,'neuropixPhase3A_kilosortChanMap.mat');

% sample rate
ops.fs = dataset.ap.meta.imSampRate;

% frequency for high pass filtering (150)
ops.fshigh = 150;

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.1;

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4];

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
ops.lam = 10;

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.95;

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/50;

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
ops.NchanTOT    = dataset.ap.meta.nSavedChans;

% common average referencing by median
ops.CAR = 1;
%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction.
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
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

fprintf('found %d good units \n', sum(rez.good>0))
%% Get sorted spikes
spike = extractrez(rez);
if ~isempty(spike)
    dataset.spike_kilosort=spike;
end

phydir = fullfile(binrootdir,[binname,'_Phy']);
if ~exist(phydir,'dir')
    mkdir(phydir);
end
rezToPhy(rez,phydir);

if exist(ops.fproc, 'file')
    delete(ops.fproc);
end
disp(['KiloSort2 Spike Sorting:    ',dataset.ap.meta.fileName,'    done.']);
%%
    function [spike]=extractrez(rez)
        spike.fs = rez.ops.fs;
        spike.time = NeuroAnalysis.Base.sample2time(rez.st3(:,1),spike.fs,dataset.secondperunit);
        spike.template = uint32(rez.st3(:,2));
        spike.amplitude = rez.st3(:,3);
        if size(rez.st3,2)>4
            spike.cluster = uint32(rez.st3(:,5));
        else
            spike.cluster = spike.template;
        end
        
        spike.ops = rez.ops;
    end
%%
    function [sspike]=splitspike(spike,range)
        sspike=spike;
        si = find(spike.time>= range(1) & spike.time< range(2));
        sspike.time=spike.time(si);
        sspike.template=spike.template(si);
        sspike.amplitude=spike.amplitude(si);
        sspike.cluster=spike.cluster(si);
    end
end