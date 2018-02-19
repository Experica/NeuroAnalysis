function [oldtest] = MergeMetadata( oldtest, newtest )
%MERGEMETADATA Merge two tests' metadata
%   Detailed explanation goes here

oldfields = fieldnames(oldtest);
newfields = fieldnames(newtest);
overwrite = intersect(oldfields, newfields);

for f = 1:length(overwrite)
    oldtest.(overwrite{f}) = newtest.(overwrite{f});
end


end
