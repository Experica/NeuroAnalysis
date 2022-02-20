function [swf] = spikefeature(waveform,fs)
%SPIKEFEATURE 1D Spike Waveform Feature
%   Detailed explanation goes here

minfs = 4.8e5; % 480kHz of 2.1us resolution
slopedt = 0.00015; % 0.15ms time interval to fit slope
if fs < minfs
    n = length(waveform);
    waveform = csaps(1:n,double(waveform),[],1:(fs/minfs):n);
    fs=minfs;
end
n = length(waveform);
slopei = 1:ceil(fs*slopedt);

[trough,troughi]=min(waveform);
[peak,peaki]=max(waveform);
amp = peak-trough;
ptr = peak/trough;
sd = (peaki-troughi)/fs;

    function [i]=searchvalueindex(start,step,stop,value)
        ssign = sign(waveform(start)-value);
        for i=start:step:stop
            if sign(waveform(i)-value)~=ssign
                break;
            end
        end
    end

% normalize to trough amplitude
ntwaveform = waveform/abs(trough);
nttrough=min(ntwaveform);
% width at half trough
lefti = searchvalueindex(troughi,-1,1,nttrough/2);
righti = searchvalueindex(troughi,1,n,nttrough/2);
lhtw = (troughi-lefti)/fs;
rhtw = (righti-troughi)/fs;
htw = (righti-lefti)/fs;
% repolarization rate
tsi = troughi+slopei;
y = ntwaveform(tsi(tsi<=n))';
X = [ones(size(y)), (1:length(y))'/fs];
c = X\y;
rprate = c(2);

% normalize to peak amplitude
npwaveform = waveform/abs(peak);
nppeak=max(npwaveform);
% width at half peak
lefti = searchvalueindex(peaki,-1,1,nppeak/2);
righti = searchvalueindex(peaki,1,n,nppeak/2);
lhpw = (peaki-lefti)/fs;
rhpw = (righti-peaki)/fs;
hpw = (righti-lefti)/fs;
% recovery rate
psi = peaki+slopei;
y = npwaveform(psi(psi<=n))';
X = [ones(size(y)), (1:length(y))'/fs];
c = X\y;
rcrate = c(2);


swf.amplitude = amp;
swf.peaktroughratio = ptr;
swf.duration = sd;
swf.ttrough = troughi/fs; % trough time

swf.lefthalftroughwidth = lhtw;
swf.righthalftroughwidth = rhtw;
swf.halftroughwidth = htw;
swf.lefthalfpeakwidth = lhpw;
swf.righthalfpeakwidth = rhpw;
swf.halfpeakwidth = hpw;
swf.repolarrate = rprate;
swf.recoverrate = rcrate;

end

