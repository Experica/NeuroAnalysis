dataexportroot='Y:/';
meta = NeuroAnalysis.IO.MetaTable(fullfile(dataexportroot,'metadata.mat'));

tests = meta.query('Subject_ID','AG2','RecordSession','V1','RecordSite','ODR1','sourceformat','SpikeGLX').Tests;
dataset = arrayfun(@(t)fullfile(dataexportroot,t.files),tests);

dataset = NeuroAnalysis.SpikeGLX.KiloSort(load(dataset{end}),'ap0');
dataset = NeuroAnalysis.SpikeGLX.KiloSort(dataset);
