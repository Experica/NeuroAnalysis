dataexportroot='Y:/';
meta = NeuroAnalysis.IO.MetaTable(fullfile(dataexportroot,'metadata.mat'));

tests = meta.query('Subject_ID','AG5','RecordSession','V1','RecordSite','1R','sourceformat','SpikeGLX').Tests;
dataset = arrayfun(@(t)fullfile(dataexportroot,t.files),tests);

% concat kilosort of all CondTests of a RecordSite
dataset = NeuroAnalysis.SpikeGLX.KiloSort(dataset);
% kilosort one of the CondTests of a RecordSite
dataset = NeuroAnalysis.SpikeGLX.KiloSort(load(dataset{end}),'ap0');

% remove high frequency drift movement
d = designfilt('lowpassiir','FilterOrder',6,'HalfPowerFrequency',0.1,'DesignMethod','butter');
imin = filtfilt(d,imin);





%% kilosort catgt flash segment for each allen probe
rootdir = 'E:\SpikeGLXData\allen-brain-observatory\visual-coding-neuropixels';
cachedir = fullfile(rootdir,'ecephys-cache');
rawdatadir = fullfile(rootdir,'raw-data');
d = 'ap0';
session_id = '839068429';
probe_id = '868929142';
probedir = fullfile(rawdatadir,session_id,probe_id);

% prepare a dataset for a probe binary file
anpdataset = struct;
anpdataset.imecindex = 0;
anpdataset.secondperunit = 1; % sec
anpdataset.filepath = fullfile(probedir,'flash.mat');
anpdataset.(d).metafile = fullfile(probedir,'flash_g0_tcat.imec0.ap.meta');
anpdataset.(d).meta = NeuroAnalysis.SpikeGLX.readmeta(anpdataset.(d).metafile);
anpdataset.(d).meta = NeuroAnalysis.SpikeGLX.parsemeta(anpdataset.(d).meta);
anpdataset.(d).meta.catgt = 1; % CatGT produced
save(anpdataset.filepath,'-struct','anpdataset','-v7.3');

anpdataset = NeuroAnalysis.SpikeGLX.KiloSort(anpdataset,d,'3');