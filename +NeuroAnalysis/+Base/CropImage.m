function [ ci ] = CropImage( img,width,height,varargin )
%CROPIMAGE Summary of this function goes here
%   Detailed explanation goes here

%% Parse arguments
p = inputParser;
addRequired(p,'img',@(x)isa(x,'numeric'));
addRequired(p,'width');
addRequired(p,'height');
addOptional(p,'position',[1,1],@(x)isa(x,'numeric'));
addOptional(p,'ispower2',true,@(x)isa(x,'logical'));
addOptional(p,'israndomposition',true,@(x)isa(x,'logical'));
parse(p,img,width,height,varargin{:});
img = p.Results.img;
width = p.Results.width;
height = p.Results.height;
position = p.Results.position;
ispower2 = p.Results.ispower2;
israndomposition = p.Results.israndomposition;
%% Check arguments
iw = size(img,2);
ih = size(img,1);
if ispower2
    width = 2^(NeuroAnalysis.Base.lastpow2(width));
    height = 2^(NeuroAnalysis.Base.lastpow2(height));
end
pmax = [ih-height+1,iw-width+1];
if any(pmax<1)
    ci=[];
    return;
end
if israndomposition
    position = [randi(pmax(1)),randi(pmax(2))];
else
    if isempty(position)
        position=[1,1];
    end
    if any(pmax-position<0)
        ci=[];
        return;
    end
end
%% Get Sub Image
ci = img(position(1):position(1)+height-1,position(2):position(2)+width-1,:);
end

