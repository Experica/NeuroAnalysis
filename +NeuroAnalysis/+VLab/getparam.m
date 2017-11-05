function [v] = getparam(paramstruct,param)
%GETPARAM Try get param value from ParamsStruct
%   Detailed explanation goes here

v=[];
names = fieldnames(paramstruct);
for i=1:length(names)
    if startsWith(names{i},param)
        v=NeuroAnalysis.VLab.tryparseparam(param, paramstruct.(names{i}));
        break;
    end
end

end

