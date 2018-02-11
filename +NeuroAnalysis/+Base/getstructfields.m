function [dst] = getstructfields(src,fields)
%GETSTRUCTFIELDS Get sub struct of fields
%   Detailed explanation goes here

dst=[];
if(isempty(src) || ~isa(src,'struct'))
    return;
end
vfn = intersect( fieldnames(src),fields);
if(isempty(vfn))
    return;
end

for i=1:length(vfn)
    dst.(vfn{i})=src.(vfn{i});
end

end

