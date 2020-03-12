function [ex] = prepare2(ex,dataset)
%PREPARE2 Prepare version 2 of Experica.Command data format
%   Detailed explanation goes here

import NeuroAnalysis.Base.* NeuroAnalysis.Experica.*
%% Pad CondTest
ctnames = fieldnames(ex.CondTest);
nct = max(cellfun(@(x)length(ex.CondTest.(x)),ctnames));
for i=1:length(ctnames)
    if length(ex.CondTest.(ctnames{i})) < nct
        ex.CondTest.(ctnames{i}){nct} = [];
    end
end
%% Parse CondTest
if isfield(ex.CondTest,'CondIndex')
    % Convert to 1-base
    ex.CondTest.CondIndex(cellfun(@isempty,ex.CondTest.CondIndex)) = {-2};
    ex.CondTest.CondIndex =cellfun(@(x)int64(x+1), ex.CondTest.CondIndex);
end
if isfield(ex.CondTest,'CondRepeat')
    ex.CondTest.CondRepeat(cellfun(@isempty,ex.CondTest.CondRepeat)) = {-1};
    ex.CondTest.CondRepeat = cellfun(@(x)int64(x), ex.CondTest.CondRepeat);
end
if isfield(ex.CondTest,'TrialIndex')
    % Convert to 1-base
    ex.CondTest.TrialIndex(cellfun(@isempty,ex.CondTest.TrialIndex)) = {-2};
    ex.CondTest.TrialIndex =cellfun(@(x)int64(x+1), ex.CondTest.TrialIndex);
end
if isfield(ex.CondTest,'TrialRepeat')
    ex.CondTest.TrialRepeat(cellfun(@isempty,ex.CondTest.TrialRepeat)) = {-1};
    ex.CondTest.TrialRepeat = cellfun(@(x)int64(x), ex.CondTest.TrialRepeat);
end
if isfield(ex.CondTest,'BlockIndex')
    % Convert to 1-base
    ex.CondTest.BlockIndex(cellfun(@isempty,ex.CondTest.BlockIndex)) = {-2};
    ex.CondTest.BlockIndex =cellfun(@(x)int64(x+1), ex.CondTest.BlockIndex);
end
if isfield(ex.CondTest,'BlockRepeat')
    ex.CondTest.BlockRepeat(cellfun(@isempty,ex.CondTest.BlockRepeat)) = {-1};
    ex.CondTest.BlockRepeat = cellfun(@(x)int64(x), ex.CondTest.BlockRepeat);
end
%% Try parse Environment Parameter
if ~isempty(ex.EnvParam)
    ex.EnvParam = tryparseparamstruct(ex.EnvParam);
end
%% Try parse Experiment Parameter
if ~isempty(ex.Param)
    ex.Param = tryparseparamstruct(ex.Param);
end
%% Init timing params, all times will be mapped to the reference timeline of DAQ Device where data are collected and time stampped
t0=54;
displaylatency = ex.Config.Display.(ex.Display_ID).Latency;
timerdriftspeed = ex.TimerDriftSpeed;
displayfallriselagdiff = ex.Config.Display.(ex.Display_ID).FallRiseLagDiff;
syncsearchradius = 20; % max jitter(ms) around perdicted sync time
%% Parse digital data
if ~isempty(dataset)
    if ~isfield(dataset,'digital')
        if isfield(dataset,'ap') && isfield(dataset.ap,'digital')
            digital = dataset.ap.digital;
            fs = dataset.ap.meta.fs;
        elseif isfield(dataset,'lf') && isfield(dataset.lf,'digital')
            digital = dataset.lf.digital;
            fs = dataset.lf.meta.fs;
        else
            digital=[];
        end
        
        if ~isempty(digital)
            dataset.digital = digital;
        end
    end
    
    if isfield(dataset,'digital')
        startsyncchidx = find(arrayfun(@(x)x.channel==ex.Config.StartSyncCh+1,dataset.digital));
        if ~isempty(startsyncchidx)
            t0=dataset.digital(startsyncchidx).time;
        end
    end
end
%%
    function searchrecover(from,to,data,latency,sr)
        names = fieldnames(ex.CondTest);
        fromnames = names(cellfun(@(x)startsWith(x,from),names));
        for c=1:length(fromnames)
            fes = ex.CondTest.(fromnames{c});
            toname= replace(fromnames{c},from,to);
            for p=1:nct
                recovered=[];
                fe = fes{p};
                if ~isempty(fe)
                    for ti=1:length(fe)
                        t=NeuroAnalysis.Experica.trysearchtime(fe(ti)+latency,data,sr);
                        recovered=[recovered,t];
                    end
                end
                if ~isempty(recovered) && all(arrayfun(@(x)isnan(x),recovered))
                    recovered=[];
                end
                ex.CondTest.(toname){p}=recovered;
            end
        end
    end
    function mergeevents(event1,event2,e2offset)
        for t=1:nct
            e1s=ex.CondTest.(event1){i};
            e2s = ex.CondTest.(event2){i}+e2offset;
            merged=e1s;
            if isempty(merged)
                if ~isempty(e2s)
                    merged=e2s;
                end
            else
                if ~isempty(e2s)
                    for ei=1:length(merged)
                        if isnan(merged(ei)) && ~isnan(e2s(ei))
                            merged(ei)=e2s(ei);
                        end
                    end
                end
            end
            if ~isempty(merged) && all(arrayfun(@(x)isnan(x),merged))
                merged=[];
            end
            ex.CondTest.(event1){i}=merged;
        end
    end
    function [uets] = uniqueeventtime(ts,es)
        if isempty(ts)
            uets= struct([]);
            return;
        end
        ues = unique(es);
        for uei=1:length(ues)
            uets.(ues{uei})=[];
        end
        for ei=1:length(es)
            uets.(es{ei})=[uets.(es{ei}),ts(ei)];
        end
    end
    function ts = firsteventtime(event)
        ts=[];
        for q=1:nct
            vs=ex.CondTest.(event){q};
            if ~isempty(vs)
                t=vs(1);
            else
                t=NaN;
            end
            ts=[ts,t];
        end
    end
%% Parse CondTest Sync Event Time, get `Command`, `Sync` and `Measure` versions of SyncEvent Times, then Combine them to get the final most accurate time if any version is corrupted.
if isfield(ex.CondTest,'Event') && isfield(ex.CondTest,'SyncEvent')
    % Parse SyncEvent Time recorded in `Command`
    synceventseq=[];synceventctidx=[];
    for i = 1:nct
        ctes = ex.CondTest.Event{i};
        ctses = ex.CondTest.SyncEvent{i};
        synceventseq=[synceventseq,ctses];
        synceventctidx=[synceventctidx,repelem(i,length(ctses))];
        usets = uniqueeventtime(time2ref(findeventtime(ctes,ctses),t0,timerdriftspeed),ctses);
        uses = fieldnames(usets);
        for j=1:length(uses)
            e = uses{j};
            ee = ['Command_',e];
            if ~isfield(ex.CondTest,ee)
                ex.CondTest.(ee)=cell(1,nct); % ensure each field has the same length
            end
            ex.CondTest.(ee){i}=usets.(e);
        end
    end
    usyncevents = unique(synceventseq);
    % Check Sync Event Data
    iseventsync = false;
    iseventmeasure = false;
    if ~isempty(dataset) && isfield(dataset,'digital')
        if ex.EventSyncProtocol.nSyncChannel==1 && ex.EventSyncProtocol.nSyncpEvent==1
            eventsyncchidx = find(arrayfun(@(x)x.channel==ex.Config.EventSyncCh+1,dataset.digital));
            eventmeasurechidx = find(arrayfun(@(x)x.channel==ex.Config.EventMeasureCh+1,dataset.digital));
            
            iseventsync = ~isempty(eventsyncchidx);
            if iseventsync
                eventsynctime = dataset.digital(eventsyncchidx).time;
                eventsyncdata = dataset.digital(eventsyncchidx).data;
                % clean noisy digital caused by logical high very close to threshold
                if length(synceventseq)~=length(eventsyncdata)
                    warning('    Clean Noisy Sync Signal.         SyncEvents/DigitalEvents:    %d/%d    ...',length(synceventseq), length(eventsyncdata));
                    if ex.PreICI ==0 && ex.SufICI==0
                        minlowdur = max(8,ex.CondDur-100);
                        minhighdur = minlowdur;
                    else
                        minlowdur = max(8,ex.PreICI+ex.SufICI-100);
                        minhighdur = max(8,ex.CondDur-100);
                    end
                    [eventsynctime,eventsyncdata] = cleannoisedigital(eventsynctime,eventsyncdata,minlowdur,minhighdur);
                end
                if length(synceventseq)==length(eventsyncdata) && all(diff(double(eventsyncdata)))
                    iseventsyncerror=false;
                else
                    iseventsyncerror=true;
                    warning('    Event Sync Error.        Events/Syncs: %d/%d,    No Flips: %d',length(synceventseq),length(eventsyncdata),find(diff(double(eventsyncdata))==0));
                end
                ex.eventsyncintegrity=~iseventsyncerror;
            end
            
            iseventmeasure = ~isempty(eventmeasurechidx);
            if iseventmeasure
                eventmeasuretime = dataset.digital(eventmeasurechidx).time;
                eventmeasuredata = dataset.digital(eventmeasurechidx).data;
                % Correct digital lag caused by asymmetrical black/white display response time,
                % here assume slow white -> black transition gives digital 0
                if displayfallriselagdiff~=0
                    fallindex = eventmeasuredata==0;
                    eventmeasuretime(fallindex) = eventmeasuretime(fallindex) - displayfallriselagdiff;
                end
                if length(synceventseq)==length(eventmeasuredata) && all(diff(double(eventmeasuredata)))
                    iseventmeasureerror=false;
                else
                    iseventmeasureerror=true;
                    warning('    Event Measure Error.        Events/Measures: %d/%d,    No Flips: %d',length(synceventseq),length(eventmeasuredata),find(diff(double(eventmeasuredata))==0));
                end
                ex.eventmeasureintegrity=~iseventmeasureerror;
            end
        end
    end
    % Parse Sync Event `Sync` Time
    if iseventsync
        if ~iseventsyncerror
            for i=1:length(synceventseq)
                se = ['Sync_',synceventseq{i}];
                if ~isfield(ex.CondTest,se)
                    ex.CondTest.(se)=cell(1,nct);
                end
                ex.CondTest.(se){synceventctidx(i)} = [ex.CondTest.(se){synceventctidx(i)},eventsynctime(i)];
            end
        else
            % Try to recover as many sync time as possible based on `Command` Time
            searchrecover('Command_','Sync_',eventsynctime,0,syncsearchradius);
        end
    end
    % Parse Sync Event `Measure` Time
    if iseventmeasure
        if ~iseventmeasureerror
            for i=1:length(synceventseq)
                me = ['Measure_',synceventseq{i}];
                if ~isfield(ex.CondTest,me)
                    ex.CondTest.(me)=cell(1,nct);
                end
                ex.CondTest.(me){synceventctidx(i)} = [ex.CondTest.(me){synceventctidx(i)},eventmeasuretime(i)];
            end
        else
            % Try to recover as many measure time as possible based on `Sync` Time
            searchrecover('Sync_','Measure_',eventmeasuretime,displaylatency,ex.Config.MaxDisplayLatencyError);
        end
    end
    
    % Try to get the most accurate and complete first Cond On/Off Time
    ismeasurecondon = isfield(ex.CondTest,'Measure_COND');
    issynccondon = isfield(ex.CondTest,'Sync_COND');
    iscommandcondon = isfield(ex.CondTest,'Command_COND');
    ismeasurecondoff = isfield(ex.CondTest,'Measure_SUFICI');
    issynccondoff = isfield(ex.CondTest,'Sync_SUFICI');
    iscommandcondoff = isfield(ex.CondTest,'Command_SUFICI');
    
    % Re-Evaluate Timing Params
    if ismeasurecondon && issynccondon
        m = firsteventtime('Measure_COND');
        s = firsteventtime('Sync_COND');
        ex.MeanSyncToMeasureDelay = nanmean(m-s);
        ex.STDSyncToMeasureDelay = nanstd(m-s);
        ex.EvalDisplayLatency = ex.MeanSyncToMeasureDelay;
        displaylatency = ex.EvalDisplayLatency;
    end
    if issynccondon && iscommandcondon
        s = firsteventtime('Sync_COND')';
        v = ref2time(firsteventtime('Command_COND')',t0,timerdriftspeed); % back to original `Command` time
        valid = ~isnan(s) & ~isnan(v);
        s = s(valid); v = v(valid);
        X = [ones(length(v),1), v];
        e = X\(s-v); % linear regression
        % intercept is the supposed fixed delay
        ex.CommandToSyncDelay = e(1);
        % slope would be the drift speed caused by different timers of Command and DAQ Device
        ex.EvalTimerDriftSpeed = e(2);
        % the fixed CommandToSyncDelay sould be very short, so the remaining would be the DAQ reference time of Command t0
        ex.t0 = ex.CommandToSyncDelay;
        % updated t0 and TimerDriftSpeed make perdicted `Sync` more accurate, so we reduce radius to get more accurate search
        syncsearchradius = 10;
    end
    
    % remap Sync Event `Command` Time based on updated timing params
    if isfield(ex,'t0')
        for i = 1:length(usyncevents)
            uce = ['Command_',usyncevents{i}];
            ex.CondTest.(uce)=cellfun(@(x)time2ref(ref2time(x,t0,timerdriftspeed),ex.t0,ex.EvalTimerDriftSpeed),ex.CondTest.(uce),'Un',0);
        end
    end
    
    % rerecover Sync Event `Sync` Time based on updated `Command` Time
    if isfield(ex,'t0') && iseventsync && iseventsyncerror
        searchrecover('Command_','Sync_',eventsynctime,0,syncsearchradius);
    end
    % rerecover Sync Event `Measure` Time based on updated DisplayLatency
    if isfield(ex,'EvalDisplayLatency') && iseventmeasure && iseventmeasureerror
        searchrecover('Sync_','Measure_',eventmeasuretime,displaylatency,ex.Config.MaxDisplayLatencyError);
    end
    
    condonversion='None';condoffversion='None';
    if ismeasurecondon
        ex.CondTest.CondOn=ex.CondTest.Measure_COND;
        condonversion='Measure';
    elseif issynccondon
        ex.CondTest.CondOn=cellfun(@(x)x+displaylatency,ex.CondTest.Sync_COND,'Un',0);
        condonversion='Sync';
    elseif iscommandcondon
        ex.CondTest.CondOn=cellfun(@(x)x+displaylatency,ex.CondTest.Command_COND,'Un',0);
    end
    if ismeasurecondoff
        ex.CondTest.CondOff=ex.CondTest.Measure_SUFICI;
        condoffversion='Measure';
    elseif issynccondoff
        ex.CondTest.CondOff=cellfun(@(x)x+displaylatency,ex.CondTest.Sync_SUFICI,'Un',0);
        condoffversion='Sync';
    elseif iscommandcondoff
        ex.CondTest.CondOff=cellfun(@(x)x+displaylatency,ex.CondTest.Command_SUFICI,'Un',0);
    end
    
    if isfield(ex.CondTest,'CondOn') && strcmp(condonversion,'Measure') && iseventmeasureerror && issynccondon
        mergeevents('CondOn','Sync_COND',displaylatency);
        condonversion='Sync';
    end
    if isfield(ex.CondTest,'CondOff') && strcmp(condoffversion,'Measure') && iseventmeasureerror && issynccondoff
        mergeevents('CondOff','Sync_SUFICI',displaylatency);
        condoffversion='Sync';
    end
    if isfield(ex.CondTest,'CondOn') && strcmp(condonversion,'Sync') && iseventsyncerror && iscommandcondon
        mergeevents('CondOn','Command_COND',displaylatency);
    end
    if isfield(ex.CondTest,'CondOff') && strcmp(condoffversion,'Sync') && iseventsyncerror && iscommandcondoff
        mergeevents('CondOff','Command_SUFICI',displaylatency);
    end
    % Try to get first Cond On/Off timing
    if isfield(ex.CondTest,'CondOn')
        ex.CondTest.CondOn=firsteventtime('CondOn');
    end
    if isfield(ex.CondTest,'CondOff')
        ex.CondTest.CondOff=firsteventtime('CondOff');
    end
    % Parse CondOff based on CondOn when PreICI == SufICI == 0
    if isfield(ex.CondTest,'CondOn') && ~isfield(ex.CondTest,'CondOff') && ex.PreICI ==0 && ex.SufICI==0
        for i=1:nct-1
            currentontime=ex.CondTest.CondOn(i);
            nextontime = ex.CondTest.CondOn(i+1);
            if (nextontime - currentontime) > (ex.CondDur+2*ex.Config.MaxDisplayLatencyError)
                ex.CondTest.CondOff(i) = currentontime + ex.CondDur;
            else
                ex.CondTest.CondOff(i)=nextontime;
            end
        end
        ex.CondTest.CondOff(nct)=ex.CondTest.CondOn(nct)+ex.CondDur;
    end
end
%% Parse CondTest Condition
if ~isempty(ex.Cond) && isfield(ex.CondTest,'CondIndex')
    % parse experiment designed condition table
    fs = fieldnames(ex.Cond);
    for i=1:length(fs)
        f=fs{i};
        ex.Cond.(f) = cellfun(@(x)tryparseparam(f,x),ex.Cond.(f),'uniformoutput',false);
    end
    
    % parse conditions been actually tested in each condtest
    for i=1:length(fs)
        ctc.(fs{i}) = cell(1, nct);
    end
    isori = ismember('Ori',fs);
    isorioffset = ismember('OriOffset',fs);
    isposition = ismember('Position',fs);
    ispositionoffset=ismember('PositionOffset',fs);
    isenvori=containsparam(ex.EnvParam,'Ori');
    isenvorioffset = containsparam(ex.EnvParam,'OriOffset');
    isenvposition=containsparam(ex.EnvParam,'Position');
    isenvpositionoffset=containsparam(ex.EnvParam,'PositionOffset');
    if containsparam(ex.EnvParam,'OriPositionOffset')
        isenvoripositionoffset=getparam(ex.EnvParam,'OriPositionOffset');
    else
        isenvoripositionoffset=false;
    end
    isorifinal=false;ispositionfinal=false;
    if (isori && (isenvorioffset || isorioffset)) || (isorioffset && (isenvori || isori))
        ctc.Ori_Final = cell(1, nct);
        isorifinal=true;
    end
    if (isposition && (isenvpositionoffset || ispositionoffset)) || (ispositionoffset && (isenvposition || isposition))
        ctc.Position_Final = cell(1, nct);
        ispositionfinal=true;
    end
    for i=1:nct
        if ex.CondTest.CondIndex(i) < 1
            continue;
        end
        
        for fi=1:length(fs)
            f=fs{fi};
            ctc.(f)(i)=ex.Cond.(f)(ex.CondTest.CondIndex(i));
        end
        % parse final orientation
        if isorifinal
            if isori
                ori = ctc.Ori{i};
            elseif isenvori
                ori = getparam(ex.EnvParam,'Ori');
            end
            if isorioffset
                orioffset = ctc.OriOffset{i};
            elseif isenvorioffset
                orioffset = getparam(ex.EnvParam,'OriOffset');
            end
            ctc.Ori_Final{i}=ori+orioffset;
        end
        % parse final position
        if ispositionfinal
            if isposition
                position=ctc.Position{i};
            elseif isenvposition
                position=getparam(ex.EnvParam,'Position');
            end
            if ispositionoffset
                positionoffset=ctc.PositionOffset{i};
            elseif isenvpositionoffset
                positionoffset = getparam(ex.EnvParam,'PositionOffset');
            end
            if isenvoripositionoffset
                theta=ori+orioffset;
                finalposition=position+positionoffset*[cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1];
            else
                finalposition=position+positionoffset;
            end
            ctc.Position_Final{i}=finalposition;
        end
    end
    ex.CondTestCond = ctc;
end
%% Standardize Experiment
ex = NeuroAnalysis.Base.StandardizeEx(ex);
end