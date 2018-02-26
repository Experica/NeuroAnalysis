function [ms] = mergestruct(varargin)
%MERGESTRUCT Merge structs with no intersected fields
%   Detailed explanation goes here

ms=[];
vsi = cellfun(@(x)~isempty(x) && isa(x,'struct') && ~isempty(fieldnames(x)),varargin);
vs = varargin(vsi);
if(isempty(vs))
    warning('No Valid Structs to Merge.');
    return;
end

fn = [];
c = [];
for i = 1:length(vs)
    fn = [fn;fieldnames(vs{i})];
    c = [c;struct2cell(vs{i})];
end

[ufn,ia,~] = unique(fn);
if length(fn) ~= length(ufn)
    warning('Fields Not Unique, use the first ones');
end

ms = cell2struct(c(ia),ufn,1);

end

