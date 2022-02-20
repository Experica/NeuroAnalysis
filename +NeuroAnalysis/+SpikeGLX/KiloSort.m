function [ dataset ] = KiloSort( dataset,d,kilosortversion )
%KILOSORT Summary of this function goes here
%   Detailed explanation goes here

if nargin==1
    d='ap';
    kilosortversion = '3';
elseif nargin==2
    kilosortversion = '3';
end

%% Concat binary files in time order into one binary file, then kilosort on it to get consistent cluster id across multiple binary files
if iscell(dataset)
    datasets = cellfun(@load,dataset,'uniformoutput',0);
    imecindex = datasets{1}.imecindex;
    ds = cellfun(@(x)['ap',x],imecindex,'uniformoutput',0);
    for k = 1:length(ds)
        d = ds{k};
        binfiles = cellfun(@(x)x.(d).meta.fileName,datasets,'uniformoutput',0);
        binfiledate = cellfun(@(x)x.(d).meta.fileDate,datasets);
        binfilensample = cellfun(@(x)x.(d).meta.nFileSamp,datasets);
        
        [~,isort]=sort(binfiledate);
        kdataset = dataset(isort);
        kdatasets = datasets(isort);
        binfiles = binfiles(isort);
        binfilensample = binfilensample(isort);
        % make concat file name
        if strcmp(kdatasets{1}.ex.sourceformat, 'Stimulator')
            binrootdir = fileparts(binfiles{1});
            concatname=cell(1,length(binfiles)+1);
            concatname{1} = kdatasets{1}.source(1:8);  % Works for AE9, not AE4 yet!
            for i=1:length(binfiles)
                concatname{i+1} = kdatasets{i}.source(10:12);
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
        concatfilepath = [concatfilepath,'.imec',imecindex{k},'.ap.bin'];
        concatmetafilepath = [concatfilepath,'.imec',imecindex{k},'.ap.meta'];
        % demux channel groups
        nchsaved = datasets{1}.(d).meta.nSavedChans;
        dmxgroup = datasets{1}.(d).meta.dmxgroup;
        if datasets{1}.(d).meta.probeversion>1
            rmdc = false;
        else
            rmdc = true;
        end
        for g = 1:length(dmxgroup)
            dmxgroup{g} = setdiff(dmxgroup{g},datasets{1}.(d).meta.excludechans);
        end
        dmxgroup(cellfun(@(x)isempty(x),dmxgroup))=[];
        % concat binary files
        if exist(concatfilepath,'file')
            disp(['Use Existing Concat Binary File:    ',concatfilepath,'    ...']);
        else
            cfid=fopen(concatfilepath,'w');
            chunksample=1e6; % 1e6 Samples = 770MB for 385Chs Int16, ≈ 34 sec for 30kS/s
            chunktestsample = round(1.5*chunksample);
            for i=1:length(binfiles)
                fprintf('Concat Binary File:    %s    ...\n',binfiles{i});
                fid=fopen(binfiles{i},'r');
                remainingsample = binfilensample(i);
                while remainingsample>0
                    if remainingsample - chunktestsample >0
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
        if exist(concatmetafilepath,'file')
            isconcatmetafile=true;
        else
            isconcatmetafile=false;
        end
        % prepare a new dataset with concat binary file
        kcdataset = struct;
        kcdataset.secondperunit = kdatasets{1}.secondperunit;
        kcdataset.filepath = kdataset;
        kcdataset.(d).meta = kdatasets{1}.(d).meta;
        kcdataset.(d).meta.fileName = concatfilepath;
        kcdataset.(d).meta.nFileSamp = sum(binfilensample);
        kcdataset.(d).meta.fromcatgt = isconcatmetafile; % CatGT produce meta file, while concat above does not
        % demuxed CAR could help to remove very fast transient noise, and then CAR in kilosort could further reduce other noise.
        kcdataset.car = 1;
        % time range [t(i), t(i+1)) for each file in the concat file
        kcdataset.binfilerange = NeuroAnalysis.Base.sample2time(cumsum([1,binfilensample]),kcdataset.(d).meta.fs,kcdataset.secondperunit);
        clear kdatasets chunkdata % reclaim memory before KiloSort
        kcdataset = NeuroAnalysis.SpikeGLX.KiloSort(kcdataset,d,kilosortversion);
        % split sorting result into each original dataset
        sf = ['spike',d(3:end),'_kilosort',kilosortversion];
        if isfield(kcdataset,sf)
            disp('Split Sorting Result Into Dataset        ...');
            for i=1:length(kdataset)
                odataset = matfile(kdataset{i},'Writable',true);
                odataset.(sf) = NeuroAnalysis.Base.splitspike(kcdataset.(sf),kcdataset.binfilerange(i:i+1));
            end
            disp('Split Sorting Result Into Dataset        Done.');
        end
    end
    return;
end

%% Choose Kilosort Version
rmpath(genpath('C:\Users\fff00\Kilosort-2.0'));
rmpath(genpath('C:\Users\fff00\Kilosort-2.5'));
rmpath(genpath('C:\Users\fff00\Kilosort'));
switch kilosortversion
    case '2'
        addpath(genpath('C:\Users\fff00\Kilosort-2.0'));
    case '25'
        addpath(genpath('C:\Users\fff00\Kilosort-2.5'));
    case '3'
        addpath(genpath('C:\Users\fff00\Kilosort'));
end

%% KiloSort ops
disp(['KiloSort ',kilosortversion,' Spike Sorting:    ',dataset.(d).meta.fileName,'    ...']);

% the binary file
ops.fbinary = dataset.(d).meta.fileName;

% if the binary file is the output of CatGT, then filtering and CAR have already been applied.
ops.fbinaryfromcatgt = false;
if isfield(dataset.(d).meta,'fromcatgt')
    ops.fbinaryfromcatgt = dataset.(d).meta.fromcatgt;
end

% the binary file folder
[binrootdir,binname,~] = fileparts(ops.fbinary);

% the probe channel map
ops.chanMap = NeuroAnalysis.SpikeGLX.neuropixelschmap(dataset.(d).meta);

% exclude channels
excludechans = dataset.(d).meta.excludechans;
if ~isempty(excludechans)
    disp(['Excluding Channels: ', num2str(excludechans)])
    ops.chanMap.connected(excludechans)=false;
end

% sample rate
ops.fs = dataset.(d).meta.fs;

% frequency for high pass filtering
ops.fshigh = 300;

% frequency for low pass filtering
% narrow spike have significent power in high frequencies(up to 10kHz), if low pass
% filtering cuts off high freq component, the spike time would be shifted
% and spike shape be widen, so here we are not doing low pass.
% ops.fslow = 7000;

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30;

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8;

% proc file on a fast SSD
ops.fproc = fullfile(binrootdir,[binname,'.kilosort',kilosortversion,'.whiten.dat']);

% time range to sort
ops.trange = [0 Inf];

% total number of channels in your recording
ops.NchanTOT = double(dataset.(d).meta.nSavedChans);

% global common average referencing by median(better to use local CAR which is not implementated yet)
ops.CAR = 1;
if isfield(dataset,'car')
    ops.CAR = dataset.car;
end

switch kilosortversion
    case '2'
        % minimum firing rate on a "good" channel (0 to skip)
        ops.minfr_goodchannels = 1/60;
        % minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
        ops.minFR = 1/60;
        % threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
        ops.Th = [12 6];
        % how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
        ops.lam = 15;
        % splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
        ops.AUCsplit = 0.9;
        % number of samples to average over (annealed from first to second value)
        ops.momentum = [20 400];
    case '3'
        % spatial smoothness constant for registration
        ops.sig = 20;
        % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
        ops.nblocks = 8;
        
        ops.Th = [10 10];
        ops.lam = 20;
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
disp(['KiloSort ',kilosortversion,' Spike Sorting:    ',dataset.(d).meta.fileName,'    done.']);

%% Get sorted spikes
disp('Extract Spike Sorting Result    ...');
switch kilosortversion
    case '2'
        spike = [];%NeuroAnalysis.SpikeGLX.extractrez2(rez,dataset.secondperunit);
    case '3'
        spike = NeuroAnalysis.SpikeGLX.extractrez3(rez,dataset.secondperunit);
end
disp('Extract Spike Sorting Result    done.');
if ~isempty(spike)
    dataset.(['spike',d(3:end),'_kilosort',kilosortversion])=spike;
end

phydir = fullfile(binrootdir,[binname,'.kilosort',kilosortversion,'.phy']);
if exist(phydir,'dir')
    rmdir(phydir,'s');
end
if ~exist(phydir,'dir')
    mkdir(phydir);
end
disp(['Save Spike Sorting Result for Phy in:    ',phydir,'    ...']);
switch kilosortversion
    case '2'
        NeuroAnalysis.SpikeGLX.rez2phy2(rez,phydir,dataset,d);
    case '3'
        NeuroAnalysis.SpikeGLX.rez2phy3(rez,phydir,dataset,d);
end
disp(['Save Spike Sorting Result for Phy in:    ',phydir,'    done.']);

end