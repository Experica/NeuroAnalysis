function [swf] = spikefeature(waveform,fs)
%SPIKEFEATURE Spike Waveform Feature
%   Detailed explanation goes here

minfs = 4.8e5; % 480kHz of 2.1us resolution
if fs < minfs
    n = length(waveform);
    waveform = csaps(1:n,waveform,[],1:(fs/minfs):n);
    fs=minfs;
end

n = length(waveform);
[trough,troughi]=min(waveform);
[peak,peaki]=max(waveform);
sd = (peaki-troughi)/fs;
ptr = peak/trough;

    function [i]=searchvalueindex(start,step,stop,value)
        ssign = sign(waveform(start)-value);
        for i=start:step:stop
            if sign(waveform(i)-value)~=ssign
                break;
            end
        end
    end

lefti = searchvalueindex(troughi,-1,1,trough/2);
righti = searchvalueindex(troughi,1,n,trough/2);
htw = (righti-lefti)/fs;
rpslope = fs*(trough/2-trough)/(righti-troughi);

lefti = searchvalueindex(peaki,-1,1,peak/2);
righti = searchvalueindex(peaki,1,n,peak/2);
hpw = (righti-lefti)/fs;
rcslope = fs*(peak/2-peak)/(righti-peaki);

swf.peaktroughratio = ptr;
swf.duration = sd;
swf.halftroughwidth = htw;
swf.halfpeakwidth = hpw;
swf.repolarslope = rpslope;
swf.recoverslope = rcslope;

end

