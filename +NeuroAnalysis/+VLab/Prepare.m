function [ ex ] = Prepare( ex, varargin )
%PREPARE Prepare VLab data
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'ex');
parse(p,ex,varargin{:});
ex = p.Results.ex;
%% Prepare data

end
