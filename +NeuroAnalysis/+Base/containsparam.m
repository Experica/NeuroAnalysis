function [r] = containsparam(paramstruct,param)
%CONTAINSPARAM Whether ParamStruct has relevant param field
%   Detailed explanation goes here

r = NeuroAnalysis.Base.matchparam(fieldnames(paramstruct),param);

end

