function [index] = MatchMetadata(test,Tests)
%MATCHMETADATA Search for matching tests
%   Detailed explanation goes here
match = false(length(Tests),1);
for t = 1:length(Tests)
    match(t) = isequal(Tests(t).filename, test.filename);
end
assert(sum(match) <= 1, 'Cannot match more than one entry')
index = find(match,1,'first');
end

