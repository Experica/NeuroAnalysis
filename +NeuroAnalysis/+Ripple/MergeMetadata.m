function [mtest] = MergeMetadata( metatable, test,index )
%MERGEMETADATA Merge test into metatable row of index
%   Detailed explanation goes here

mtest = metatable.Tests(index);

oldfields = fieldnames(metatable.Tests);
newfields = fieldnames(test);
overwrite = intersect(oldfields, newfields);

for i = 1:length(overwrite)
    mtest.(overwrite{i}) = test.(overwrite{i});
end

end
