function [mtest] = MergeMetadata( metatable, test,index )
%MERGEMETADATA Merge test into metatable row of index
%   Detailed explanation goes here

mtest = metatable.Tests(index);
fields = fieldnames(test);

isoimerge = strcmp(mtest.sourceformat,'OI');
if isoimerge
    fields(cellfun(@(x)strcmp(x,'files'),fields))=[];
    fields(cellfun(@(x)strcmp(x,'sourceformat'),fields))=[];
end

for i = 1:length(fields)
    mtest.(fields{i}) = test.(fields{i});
end

if isoimerge && ~ismember(test.files,mtest.files)
    mtest.files=[mtest.files;test.files];
end

end