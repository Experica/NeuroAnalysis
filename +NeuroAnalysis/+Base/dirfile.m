function [ filedir,filename,filedate ] = dirfile( rootdir,nametemplate,sortmethod )
%DIRFILE Get file directory, name and date under a root directory
%   Detailed explanation goes here

if nargin==2
    sortmethod='descend';
elseif nargin==1
    nametemplate='*';
    sortmethod='descend';
end

fs = dir(fullfile(rootdir,nametemplate));
fsname = arrayfun(@(i)i.name,fs,'UniformOutput',false);
fsdir = arrayfun(@(i)i.folder,fs,'UniformOutput',false);
fsdatenum = arrayfun(@(i)i.datenum,fs);
[~,sorti] = sort(fsdatenum,sortmethod);

filename = fsname(sorti);
filedir = fsdir(sorti);
filedate = fsdatenum(sorti);

end

