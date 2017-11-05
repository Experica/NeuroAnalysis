function [r] = containsparam(paramstruct,param)
%CONTAINSPARAM Whether ParamsStruct has param field
%   Detailed explanation goes here

r=any(startsWith(fieldnames(paramstruct),param));

end

