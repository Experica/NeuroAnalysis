function [v] = getparam(paramstruct,param,ignorecase)
%GETPARAM Try get param value from ParamStruct
%   Detailed explanation goes here

if nargin < 3
    ignorecase = false;
end

import NeuroAnalysis.Base.matchparam

v=[];
names = fieldnames(paramstruct);
for i=1:length(names)
    if matchparam(names{i}, param, ignorecase)
        v=paramstruct.(names{i});
        break;
    end
end

end

