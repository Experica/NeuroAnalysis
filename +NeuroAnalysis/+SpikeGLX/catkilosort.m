dataexportroot='Y:/';
meta = NeuroAnalysis.IO.MetaTable(fullfile(dataexportroot,'metadata.mat'));

tests = meta.query('Subject_ID','AG1','RecordSession','V1','RecordSite','ODL1','sourceformat','SpikeGLX').Tests;
dataset = arrayfun(@(t)fullfile(dataexportroot,t.files),tests);

% concat kilosort of all CondTests of a RecordSite
dataset = NeuroAnalysis.SpikeGLX.KiloSort(dataset);
% kilosort one of the CondTests of a RecordSite
dataset = NeuroAnalysis.SpikeGLX.KiloSort(load(dataset{end}),'ap0');

% remove high frequency drift movement
d = designfilt('lowpassiir','FilterOrder',6,'HalfPowerFrequency',0.1,'DesignMethod','butter');
imin = filtfilt(d,imin);





% prepare a new dataset for a Allen Neuropixels probe binary file
anpdataset = struct;
d = 'ap0';
anpdataset.secondperunit = 0.001; % ms
anpdataset.filepath = '';
anpdataset.(d).meta = NeuroAnalysis.SpikeGLX.readmeta('X:\allen-brain-observatory\visual-coding-neuropixels\raw-data\839068429\868929138\flash_g0_t0.imec0.ap.meta');
anpdataset.(d).meta = NeuroAnalysis.SpikeGLX.parsemeta(anpdataset.(d).meta);
anpdataset.(d).meta.catgt = 0; % CatGT produce meta file, while concat above does not
% demuxed CAR could help to remove very fast transient noise, and then CAR in kilosort could further reduce other noise.
anpdataset.car = 1;

anpdataset = NeuroAnalysis.SpikeGLX.KiloSort(anpdataset,d,'3');