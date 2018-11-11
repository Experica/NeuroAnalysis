function [dataset] = loaddataset(filepath,vars)
%LOADDATASET Summary of this function goes here
%   Detailed explanation goes here

dataset=[];
finfo=whos('-file',filepath);
varnames = {finfo.name};
vartypes = {finfo.class};
structindex=arrayfun(@(t)strcmp(t{:},'struct'),vartypes);
structvarnames=varnames(structindex);
nonestructvarnames=varnames(~structindex);
basevars={'ex','imagehead'};

if nargin==1
    vars=varnames;
elseif isa(vars,'cell')
    vars=union(nonestructvarnames, intersect(structvarnames,[vars,basevars]));
else
    warning('Invalid vars, should be cell of string.');
    return;
end

dataset=load(filepath,vars{:});
end