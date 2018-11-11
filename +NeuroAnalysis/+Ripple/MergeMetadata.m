function [mtest] = MergeMetadata( metatable, test,index )
%MERGEMETADATA Merge test into metatable row of index
%   Detailed explanation goes here

mtest = metatable.Tests(index);
fields = fieldnames(test);

for i = 1:length(fields)
    mtest.(fields{i}) = test.(fields{i});
end

end