function [ result ] = MetadataExtract( dataset,exportpath,sourceformat,isparallel,varargin )
%METADATAUPDATE Update metadata for exported files
%   Detailed explanation goes here

%% Batch update
if isa(dataset,'cell')
    funlist=repelem({'NeuroAnalysis.IO.MetadataExtract'},length(dataset));
    vararginlist = arrayfun(@(i)[i,{exportpath,sourceformat,isparallel},varargin],dataset,'UniformOutput',false);
    result = NeuroAnalysis.Base.ApplyFunctions(funlist,vararginlist,isparallel);
    return;
end

%% Open dataset and read metadata
result.status = false;
result.source = '';
result.meta = struct([]);
try
    if isa(dataset,'char')
        disp(['Reading Dataset:    ',dataset,'    ...']);
        result.source = dataset;
        if ~strcmp(sourceformat,'Unknown')
            filename = NeuroAnalysis.Base.filenamenodirext(dataset);
            filepath = fullfile(exportpath,[filename,'.mat']);
            datafile = load(filepath); 
            dataset = datafile.dataset;
        else
            [~,filename,ext]=fileparts(dataset);
            filepath=fullfile(exportpath,[filename,ext]);
        end
        disp('Reading Dataset:    Done.');
    else
        filepath = dataset.filepath;
    end
    % Copy metadata for this test
    test = struct;
    test.filepath = filepath;
    test.sourceformat = sourceformat;
    test.dateadded = now;
    if ~isempty(dataset) && isfield(dataset, 'ex')
        fields = {{'ID','ID',@string}, ...
            {'Subject_ID','Subject_ID',@string}, ...
            {'File_ID','File_ID',@double}, ...
            {'RecordSite','RecordSite',@string}, ...
            {'RecordSession','RecordSession',@string}};
        test = NeuroAnalysis.Base.copyStructFields(dataset.ex, test, fields);
    end
    result.status = true;
    result.source = filepath;
    result.meta = test;
catch ME
    result.message = ME.message;
    warning([ME.identifier,': ',ME.message]);
    for i=1:length(ME.stack)
        ME.stack(i);
    end
end
end
