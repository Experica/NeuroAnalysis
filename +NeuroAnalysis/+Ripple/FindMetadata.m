function [i] = FindMetadata(metatable,test,range)
%FINDMETADATA Search for matching tests in range, return indices
%   Detailed explanation goes here

for i = range
    if strcmp(metatable.Tests(i).filename, test.filename)
        return;
    end
end

i = [];

end

