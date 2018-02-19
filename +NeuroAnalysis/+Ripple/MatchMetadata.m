function [index] = MatchMetadata(test,Tests)
%MATCHMETADATA Search for matching tests
%   Detailed explanation goes here
for index = 1:length(Tests)
    if strcmp(Tests(index).filename, test.filename)
        return;
    end
end
index = [];
return
end

