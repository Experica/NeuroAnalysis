function [chmap] = neuropixelschmap(meta)
%NEUROPIXELSCHMAP generate all connected channels map
%   Detailed explanation goes here

if meta.probeversion <= 1
    chmap.name = 'Neuropixels Phase3A/3B/1.0';
    chmap.connected = true(size(meta.roch));
    chmap.shankind = meta.savedshanks;
    chmap.chanMap0ind = double(meta.robank*meta.acqApLfSy(1) + meta.roch);
    chmap.chanMap = chmap.chanMap0ind+1;
    dx = meta.probespacing(1);dy=meta.probespacing(2);
    cols = double(meta.savedcols);rows = double(meta.savedrows);
    
    % checkboard
    chmap.xcoords = arrayfun(@(r,c)(dx/2 + (c-1)*dx) - abs(mod(r,2)-1)*dx/2,rows,cols);
    chmap.ycoords = (rows-1)*dy;
    % regular
    %             chmap.xcoords = (cols-1)*dx;
    %             chmap.ycoords = (rows-1)*dy;
end

end

