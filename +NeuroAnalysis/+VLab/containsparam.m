function [r] = containsparam(paramstruct,param)
%CONTAINSPARAM Whether ParamStruct has relevant param field
%   Detailed explanation goes here

r=any(startsWith(fieldnames(paramstruct),param));

end

