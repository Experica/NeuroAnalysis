function [swf] = spikefeature(waveform,fs)
%SPIKEFEATURE 1D Spike Waveform Feature
%   Detailed explanation goes here

minfs = 4.8e5; % 480kHz of 2.1us resolution
slopedt = 0.0001; % 0.1ms time interval to fit slope
if fs < minfs
    n = length(waveform);
    waveform = csaps(1:n,double(waveform),[],1:(fs/minfs):n);
    fs=minfs;
end
n = length(waveform);
slopei = 1:ceil(fs*slopedt);

% normalize to trough amplitude
waveform = waveform/abs(min(waveform));
[trough,troughi]=min(waveform);
[peak,peaki]=max(waveform);
sd = (peaki-troughi)/fs;
ptr = peak/trough;
amp = peak-trough;

    function [i]=searchvalueindex(start,step,stop,value)
        ssign = sign(waveform(start)-value);
        for i=start:step:stop
            if sign(waveform(i)-value)~=ssign
                break;
            end
        end
    end

% width at half trough
lefti = searchvalueindex(troughi,-1,1,trough/2);
righti = searchvalueindex(troughi,1,n,trough/2);
htw = (righti-lefti)/fs;
% repolarization rate
tsi = troughi+slopei;
y = waveform(tsi(tsi<=n))';
X = [ones(size(y)), (1:length(y))'/fs];
c = X\y;
rprate = c(2);

% width at half peak
lefti = searchvalueindex(peaki,-1,1,peak/2);
righti = searchvalueindex(peaki,1,n,peak/2);
hpw = (righti-lefti)/fs;
% recovery rate
psi = peaki+slopei;
y = waveform(psi(psi<=n))';
X = [ones(size(y)), (1:length(y))'/fs];
c = X\y;
rcrate = c(2);


swf.ttrough = troughi/fs; % trough time
swf.peaktroughratio = ptr;
swf.duration = sd;
swf.amplitude = amp;
swf.halftroughwidth = htw;
swf.halfpeakwidth = hpw;
swf.repolarrate = rprate;
swf.recoverrate = rcrate;

end

