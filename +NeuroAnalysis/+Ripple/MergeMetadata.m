function [merged] = MergeMetadata( test1, test2 )
%MERGEMETADATA Merge two tests' metadata
%   Detailed explanation goes here

t1fields = fieldnames(test1);
t2fields = fieldnames(test2);
t1missing = setdiff(t2fields, t1fields);
t2missing = setdiff(t1fields, t2fields);
if test1.dateadded > test2.dateadded
    merged = test1;
    for f = 1:length(t1missing)
        merged.(t1missing{f}) = test2.(t1missing{f});
    end
else
    merged = test2;
    for f = 1:length(t2missing)
        merged.(t2missing{f}) = test1.(t2missing{f});
    end
end

end
