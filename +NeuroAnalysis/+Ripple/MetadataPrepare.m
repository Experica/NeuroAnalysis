function [ result ] = PrepareMetadata( dataset,callbackresult )
%PREPAREMETADATA Update metadata for exported files
%   Detailed explanation goes here

%% Read metadata
result.status = false;
result.source = '';
result.meta = struct([]);
test = struct;
test.filepath = dataset.filepath;
[test.filedir, test.filename, test.ext] = fileparts(test.filepath);
test.sourceformat = dataset.sourceformat;
test.dateadded = now;
if ~isempty(dataset) && isfield(dataset, 'ex')
    fields = {'ID', 'Subject_ID', 'File_ID',...
        'RecordSite','RecordSession'};
    test = NeuroAnalysis.Base.copyStructFields(dataset.ex, test, fields);
end
if ~isempty(callbackresult) && isa(callbackresult,'struct')
    test = NeuroAnalysis.Base.copyStructFields(callbackresult, test, ...
        fieldnames(callbackresult));
end
test.key = 'filename';
result.status = true;
result.source = dataset.filepath;
result.meta = test;
end
