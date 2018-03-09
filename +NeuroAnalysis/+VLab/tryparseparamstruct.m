function [paramstruct] = tryparseparamstruct(paramstruct)
%TRYPARSEPARAMSTRUCT Summary of this function goes here
%   Detailed explanation goes here

fn = fieldnames(paramstruct);
oldfieldname={};
for i=1:length(fn)
    p = fn{i};
    paramstruct.(p) =NeuroAnalysis.VLab.tryparseparam(p,paramstruct.(p));
    if contains(p,'0x40')
        paramstruct.(replace(p,'0x40','@')) = paramstruct.(p);
        oldfieldname=[oldfieldname,{p}];
    end
end
paramstruct=rmfield(paramstruct,oldfieldname);

end

