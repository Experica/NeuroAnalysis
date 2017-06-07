function [ names ] = filenamenodirext( files )
%FILENAMENODIREXT Get file name without file directory and extension
%   Detailed explanation goes here

if isa(files,'cell')
    n = length(files);
    names = cell(n,1);
    for i=1:n
        [~,name,~]=fileparts(files{i});
        names{i}=name;
    end
elseif isa(files,'char')
    [~,names,~]=fileparts(files);
else
    names=[];
end

end

