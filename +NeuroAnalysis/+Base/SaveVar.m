function [  ] = SaveVar( pathitem,var,varname )
%SAVEVAR Save a varable to a MAT-File located in the same folder of
%pathitem.
%   Detailed explanation goes here

eval([varname '=var;']);
path = fileparts(which(pathitem));
filepath = fullfile(path,[varname,'.mat']);
save(filepath,varname,'-v7.3');

end

