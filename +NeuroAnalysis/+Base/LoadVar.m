function [ var ] = LoadVar( pathitem, filename, varname )
%LOADVAR Load a varable from a MAT-File located in the same folder of
%pathitem
%   Detailed explanation goes here

if nargin==2
    [~, varname, ~] = fileparts(filename);
end

path = fileparts(which(pathitem));
filepath = fullfile(path,[filename,'.mat']);

if exist(filepath,'file')==2
    var = load(filepath,varname);
    var = var.(varname);
else
    var=[];
end

end

