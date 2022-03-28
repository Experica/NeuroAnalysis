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