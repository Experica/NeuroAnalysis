function ex = adjustStimTimes( ex, dataset, visstimconfig )
%AdjustStimTimes Adjust stim times to pick the most accurate ones

stimTimesPTB = [];
stimTimesPhotodiode = []; 
stimTimesParallel = []; 
ParallelDCh = visstimconfig.ParallelDCh;
PhotodiodeDCh = visstimconfig.PhotodiodeDCh;
NMarkPerCond = visstimconfig.NMarkPerCond;
SearchRadius = visstimconfig.MarkSearchRadius;
TimerDriftError = visstimconfig.TimerDriftError;



% Collect digital marks if they exist
nct = length(ex.CondTest.CondIndex);
photodiodeError=true;
parallelError=true;
if ~isempty(dataset) && isfield(dataset,'digital')
    parallelIdx = find(arrayfun(@(x)x.channel==ParallelDCh,dataset.digital));
    if ~isempty(parallelIdx)
        stimTimesParallel = dataset.digital(parallelIdx).time';
        stimValuesParallel = dataset.digital(parallelIdx).data';
        if NMarkPerCond*nct == length(stimTimesParallel) && all(diff(stimValuesParallel))
            parallelError = false;
        end
    end
    photoIdx = find(arrayfun(@(x)x.channel==PhotodiodeDCh,dataset.digital));
    if ~isempty(photoIdx)
        stimTimesPhotodiode = dataset.digital(photoIdx).time';
        stimValuesPhotodiode = dataset.digital(photoIdx).data';
        if NMarkPerCond*nct == length(stimTimesPhotodiode) && all(diff(stimValuesPhotodiode))
            photodiodeError = false;
        end
    end
end

% Collect PTB/matlab times from ex structure
if isfield(ex.CondTest,'StimOn') && isfield(ex.CondTest,'StimOff')
    nMarks = length(ex.CondTest.StimOn)*2;
    stimTimesPTB = nan(nMarks, 1);
    stimTimesPTB(1:2:nMarks-1) = ex.CondTest.StimOn;
    stimTimesPTB(2:2:nMarks) = ex.CondTest.StimOff;
    
    % Correct for timer drift (roughly) and ripple start offset
    X = [ones(nMarks,1) stimTimesPTB];
    stimTimesPTB = stimTimesPTB + X*TimerDriftError + ex.t0;
end
    
% Fix if broken
rippleFirst = min(min([stimTimesParallel; Inf]), min([stimTimesPhotodiode; Inf]));
if rippleFirst < Inf && rippleFirst > min(stimTimesPTB(:,1)) + 1
    % Not lined up within 1 second usually means the start time is wrong
    stimTimesPTB = stimTimesPTB + rippleFirst - min(stimTimesPTB);
end

% If photodiode is intact, use it
if ~photodiodeError
    
    corrected = stimTimesPhotodiode;
    on = corrected(1:2:end);
    off = corrected(2:2:end);
    
% If not, correct them
elseif ~isempty(stimTimesPhotodiode)
    
    corrected = correctTimes(stimTimesPTB+ex.Latency, ...
        stimTimesPhotodiode, stimValuesPhotodiode, SearchRadius);
    
    errors = ismissing(corrected);    
    corrected(errors) = stimTimesPTB(errors) + ex.Latency;
    on = corrected(1:2:end);
    off = corrected(2:2:end);

% If no photodiode, use parallel
elseif ~parallelError
    
    corrected = stimTimesParallel + ex.Latency;
    on = corrected(1:2:end);
    off = corrected(2:2:end);
    
% Correct parallel if broken
elseif ~isempty(stimTimesParallel)
    corrected = correctTimes(stimTimesPTB + ex.Latency, ...
        stimTimesParallel + ex.Latency, stimValuesParallel, SearchRadius);
    
    errors = ismissing(corrected);    
    corrected(errors) = stimTimesPTB(errors) + ex.Latency;
    on = corrected(1:2:end);
    off = corrected(2:2:end);
    
% No corrections, don't have enough data    
else
    
    on = stimTimesPTB(1:2:end) + ex.Latency;
    off = stimTimesPTB(2:2:end) + ex.Latency;
    
end
    

% Correct experiments with short pulse stim marks (should be none for
% VisStim experiments)
if min(off - on) < SearchRadius % Search radius sets the pulse minimum
    off = on + stimTimesPTB(2:2:end) - stimTimesPTB(1:2:end);
end

% None-ICI Mark Mode
if ex.PreICI==0 && ex.SufICI==0 && ~isempty(on) && ~isempty(off)
    for i=1:nct-1
        currentontime=on(i);
        nextontime = on(i+1);
        if (nextontime - currentontime) > (ex.CondDur+2*SearchRadius)
            off(i) = currentontime + ex.CondDur;
        else
            off(i)=nextontime;
        end
    end
    off(end)=on(end)+ex.CondDur;
end

difference = on - stimTimesPTB(1:2:end);
latency = mean(difference);
variation = std(difference);

ex.CondTest.CondOn = on;
ex.CondTest.CondOff = off;
ex.latency = latency;
ex.variation = variation;

end


function corrections = correctTimes(reference, newTimes, ...
    newValues, searchRadius)

corrections = nan(length(reference),1);
i = 1; % next new timestamp

% Correct each timestamp t in the reference series
for t = 1:size(reference,1)
    if i > size(newTimes,1)
        break;
    end
    j = 0;
    found = 0;
    % Search starting at i until t > reference(t) + search radius
    while ~found && i+j <= size(newTimes,1) && ...
            newTimes(i+j) < reference(t) + searchRadius
        if  newTimes(i+j) > reference(t) - searchRadius && ...
                mod(t,2) == ~~newValues(i+j)
            corrections(t) = newTimes(i+j);
            found = 1;
            break;
        end
        j = j + 1;
    end
    if found || (i+j <= size(newTimes, 1) && ...
            newTimes(i+j) >= reference(t) + searchRadius)
        i = i + j; % no need to search any t < reference(t)
    end
end

end
