function [ result ] = PrepareMetadata( dataset,callbackresult )
%PREPAREMETADATA Update metadata for exported files
%   Detailed explanation goes here

%% Read metadata
result.status = false;
result.source = '';
result.meta = struct([]);
test = struct;
% Required fields
test.files = {dataset.filepath};
test.sourceformat = dataset.sourceformat;
test.dateadded = now;
% Ripple-specific fields
[~, test.filename, ~] = fileparts(dataset.filepath);
if ~isempty(dataset) && isfield(dataset, 'ex')
    fields = {'ID', 'Subject_ID', 'File_ID',...
        'RecordSite','RecordSession'};
    test = NeuroAnalysis.Base.copyStructFields(dataset.ex, test, fields, @(x)char(string(x)));
end
if ~isempty(callbackresult) && isa(callbackresult,'struct')
    test = NeuroAnalysis.Base.copyStructFields(callbackresult, test, ...
        fieldnames(callbackresult));
end
result.status = true;
result.source = dataset.filepath;
result.meta = test;
end

