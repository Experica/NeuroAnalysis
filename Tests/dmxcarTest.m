import NeuroAnalysis.Base.dmxcar
%% simulate data
nch = 385;
nsample = 1000;
nchmx = 12;
ndmx = floor(nch/nchmx);

data = randn(nch,nsample);
data = data + 10*randn(nch,1);
for g = 0:ndmx-1
    for i = 1:nchmx
        ch = g*nchmx + i;
        for j = 1:nsample
            data(ch,j) = data(ch,j) + 10*sin(2*pi*(0.01*j-0.1*i));
        end
    end
end
data = int16(data);

pn=2*nchmx;
figure
for i =1:pn
    subplot(pn,1,i);
    plot(data(i,:))
    axis off
end

figure
imagesc(data)
%% Ordinary CAR
odata = data - mean(data,2,'native');
odata = odata - median(odata,1);

figure
for i =1:pn
    subplot(pn,1,i);
    plot(odata(i,:))
    axis off
end

figure
imagesc(odata)

%% Global Demuxed CAR
vch = 1:nch;
dmxgroup = cell(nchmx,1);
for g = 1:nchmx
    dmxgroup{g} = intersect(g + (0:nchmx:nch-nchmx),vch);
end
dmxgroup(cellfun(@(x)isempty(x),dmxgroup))=[];

ddata = dmxcar(data,dmxgroup);

figure
for i =1:pn
    subplot(pn,1,i);
    plot(ddata(i,:))
    axis off
end

figure
imagesc(ddata)
%% performance test
tic
ddata = dmxcar(data,dmxgroup);
toc