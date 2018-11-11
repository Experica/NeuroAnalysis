function [i] = FindMetadata(metatable,test,range)
%FINDMETADATA Search for matching tests in range, return indices
%   Detailed explanation goes here

i = metatable.iquery({'sourceformat'}, {'OI'},range);

end