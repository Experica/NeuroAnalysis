function [paramstruct] = tryparseparamstruct(paramstruct)
%TRYPARSEPARAMSTRUCT Summary of this function goes here
%   Detailed explanation goes here

fn = fieldnames(paramstruct);
for i=1:length(fn)
    p = fn{i};
    paramstruct.(p) =NeuroAnalysis.Experica.tryparseparam(p,paramstruct.(p));
end

end

