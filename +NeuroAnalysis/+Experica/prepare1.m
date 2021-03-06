function [ex] = prepare1(ex,dataset)
%PREPARE1 Prepare version 1 of Experica.Command data format
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
%% Parse t0
ex.t0=0;
if ~isempty(dataset) && isfield(dataset,'digital')
    startsyncchidx = find(arrayfun(@(x)x.channel==ex.Config.StartSyncCh+1,dataset.digital));
    if ~isempty(startsyncchidx)
        ex.t0=dataset.digital(startsyncchidx).time;
    end
end
%% Parse CondTest
if isfield(ex.CondTest,'CondIndex')
    % Convert to 1-base
    ex.CondTest.CondIndex(cellfun(@isempty,ex.CondTest.CondIndex)) = {-2};
    ex.CondTest.CondIndex =cellfun(@(x)int32(x+1), ex.CondTest.CondIndex);
end
if isfield(ex.CondTest,'CondRepeat')
    ex.CondTest.CondRepeat(cellfun(@isempty,ex.CondTest.CondRepeat)) = {-1};
    ex.CondTest.CondRepeat = cellfun(@(x)int32(x), ex.CondTest.CondRepeat);
end
if isfield(ex.CondTest,'TrialIndex')
    % Convert to 1-base
    ex.CondTest.TrialIndex(cellfun(@isempty,ex.CondTest.TrialIndex)) = {-2};
    ex.CondTest.TrialIndex =cellfun(@(x)int32(x+1), ex.CondTest.TrialIndex);
end
if isfield(ex.CondTest,'TrialRepeat')
    ex.CondTest.TrialRepeat(cellfun(@isempty,ex.CondTest.TrialRepeat)) = {-1};
    ex.CondTest.TrialRepeat = cellfun(@(x)int32(x), ex.CondTest.TrialRepeat);
end
if isfield(ex.CondTest,'BlockIndex')
    % Convert to 1-base
    ex.CondTest.BlockIndex(cellfun(@isempty,ex.CondTest.BlockIndex)) = {-2};
    ex.CondTest.BlockIndex =cellfun(@(x)int32(x+1), ex.CondTest.BlockIndex);
end
if isfield(ex.CondTest,'BlockRepeat')
    ex.CondTest.BlockRepeat(cellfun(@isempty,ex.CondTest.BlockRepeat)) = {-1};
    ex.CondTest.BlockRepeat = cellfun(@(x)int32(x), ex.CondTest.BlockRepeat);
end
%% Parse CondTest Sync Event Timing
    function searchrecover(from,to,data,latency,sr)
        names = fieldnames(ex.CondTest);
        fromnames = names(cellfun(@(x)startsWith(x,from),names));
        for c=1:length(fromnames)
            fct = ex.CondTest.(fromnames{c});
            te= replace(fromnames{c},from,to);
            for p=1:nct
                recovered=[];
                ctfets = fct{p};
                if ~isempty(ctfets)
                    for ti=1:length(ctfets)
                        t=NeuroAnalysis.Experica.trysearchtime(ctfets(ti)+latency,data,sr);
                        recovered=[recovered,t];
                    end
                end
                if ~isempty(recovered) && all(arrayfun(@(x)isnan(x),recovered))
                    recovered=[];
                end
                ex.CondTest.(te){p}=recovered;
            end
        end
    end
    function combinetiming(event,combineevent,latency)
        for i=1:nct
            ctts=ex.CondTest.(event){i};
            cttcs = ex.CondTest.(combineevent){i}+latency;
            if isempty(ctts)
                combined=[];
                if ~isempty(cttcs)
                    combined=cttcs;
                end
            else
                combined=ctts;
                if ~isempty(cttcs)
                    for j=1:length(ctts)
                        if isnan(ctts(j)) && ~isnan(cttcs(j))
                            combined(j)=cttcs(j);
                        end
                    end
                end
            end
            if ~isempty(combined) && all(arrayfun(@(x)isnan(x),combined))
                combined=[];
            end
            ex.CondTest.(event){i}=combined;
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
if isfield(ex.CondTest,'Event') && isfield(ex.CondTest,'SyncEvent')
    % Parse Sync Event Experica.Command Timing
    ses=[];sectidx=[];
    for i = 1:nct
        ctes = ex.CondTest.Event{i};
        ctses = ex.CondTest.SyncEvent{i};
        ses=[ses,ctses];
        sectidx=[sectidx,repelem(i,length(ctses))];
        usets = uniqueeventtime(arrayfun(@(x)toreftime(x,ex.t0,ex.TimerDriftSpeed), findeventtime(ctes,ctses)),ctses);
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
    % Check Sync Event Data
    isdineventsync = false;
    isdineventmeasure = false;
    if ~isempty(dataset) && isfield(dataset,'digital')
        if ex.EventSyncProtocol.nSyncChannel==1 && ex.EventSyncProtocol.nSyncpEvent==1
            eventsyncchidx = find(arrayfun(@(x)x.channel==ex.Config.EventSyncCh+1,dataset.digital));
            eventmeasurechidx = find(arrayfun(@(x)x.channel==ex.Config.EventMeasureCh+1,dataset.digital));
            
            isdineventsync = ~isempty(eventsyncchidx);
            if isdineventsync
                dineventsynctime = dataset.digital(eventsyncchidx).time;
                dineventsyncdata = dataset.digital(eventsyncchidx).data;
                if length(ses)==length(dineventsyncdata) && all(diff(dineventsyncdata))
                    isdineventsyncerror=false;
                else
                    isdineventsyncerror=true;
                end
                ex.eventsyncintegrity=~isdineventsyncerror;
            end
            
            isdineventmeasure = ~isempty(eventmeasurechidx);
            if isdineventmeasure
                dineventmeasuretime = dataset.digital(eventmeasurechidx).time;
                dineventmeasuredata = dataset.digital(eventmeasurechidx).data;
                if length(ses)==length(dineventmeasuredata) && all(diff(dineventmeasuredata))
                    isdineventmeasureerror=false;
                else
                    isdineventmeasureerror=true;
                end
                ex.eventmeasureintegrity=~isdineventmeasureerror;
            end
        end
    end
    % Parse Sync Event Timing
    if isdineventsync
        if ~isdineventsyncerror
            for i=1:length(ses)
                se = ['Sync_',ses{i}];
                if ~isfield(ex.CondTest,se)
                    ex.CondTest.(se)=cell(1,nct);
                end
                ex.CondTest.(se){sectidx(i)} = [ex.CondTest.(se){sectidx(i)},dineventsynctime(i)];
            end
        else
            % Try to recover as many sync timing as possible based on Experica.Command Timing
            searchrecover('Command_','Sync_',dineventsynctime,0,ex.Config.MaxDisplayLatencyError);
        end
    end
    % Parse Sync Event Measure Timing
    if isdineventmeasure
        if ~isdineventmeasureerror
            for i=1:length(ses)
                me = ['Measure_',ses{i}];
                if ~isfield(ex.CondTest,me)
                    ex.CondTest.(me)=cell(1,nct);
                end
                ex.CondTest.(me){sectidx(i)} = [ex.CondTest.(me){sectidx(i)},dineventmeasuretime(i)];
            end
        else
            % Try to recover as many measure timing as possible based on Sync Timing
            searchrecover('Sync_','Measure_',dineventmeasuretime,ex.DisplayLatency,ex.Config.MaxDisplayLatencyError);
        end
    end
    % Try to get the most accurate and complete Cond On/Off Time
    ismeasurecondon = isfield(ex.CondTest,'Measure_COND');
    issynccondon = isfield(ex.CondTest,'Sync_COND');
    iscommandcondon = isfield(ex.CondTest,'Command_COND');
    ismeasurecondoff = isfield(ex.CondTest,'Measure_SUFICI');
    issynccondoff = isfield(ex.CondTest,'Sync_SUFICI');
    iscommandcondoff = isfield(ex.CondTest,'Command_SUFICI');
    
    % Calculate measured display latency and timer drift speed
    if ismeasurecondon && issynccondon
        m = firsteventtime('Measure_COND');
        s = firsteventtime('Sync_COND');
        ex.MeasureMeanDisplayLatency = nanmean(m-s);
        ex.MeasureSTDDisplayLatency = nanstd(m-s);
    end
    if issynccondon && iscommandcondon
        s = firsteventtime('Sync_COND')';
        v = arrayfun(@(x)reftotime(x,ex.t0,ex.TimerDriftSpeed),firsteventtime('Command_COND')');
        valid = ~isnan(s) & ~isnan(v);
        s = s(valid); v = v(valid);
        X = [ones(length(s),1), v];
        e = X\(s-v);
        ex.MeasureTimerDriftSpeed = e(2);
    end
    
    condonversion='None';condoffversion='None';
    if ismeasurecondon
        ex.CondTest.CondOn=ex.CondTest.Measure_COND;
        condonversion='Measure';
    elseif issynccondon
        ex.CondTest.CondOn=cellfun(@(x)x+ex.DisplayLatency,ex.CondTest.Sync_COND,'Un',0);
        condonversion='Sync';
    elseif iscommandcondon
        ex.CondTest.CondOn=cellfun(@(x)x+ex.DisplayLatency,ex.CondTest.Command_COND,'Un',0);
    end
    if ismeasurecondoff
        ex.CondTest.CondOff=ex.CondTest.Measure_SUFICI;
        condoffversion='Measure';
    elseif issynccondoff
        ex.CondTest.CondOff=cellfun(@(x)x+ex.DisplayLatency,ex.CondTest.Sync_SUFICI,'Un',0);
        condoffversion='Sync';
    elseif iscommandcondoff
        ex.CondTest.CondOff=cellfun(@(x)x+ex.DisplayLatency,ex.CondTest.Command_SUFICI,'Un',0);
    end
    
    if isfield(ex.CondTest,'CondOn') && strcmp(condonversion,'Measure') && isdineventmeasureerror && issynccondon
        combinetiming('CondOn','Sync_COND',ex.DisplayLatency);
        condonversion='Sync';
    end
    if isfield(ex.CondTest,'CondOff') && strcmp(condoffversion,'Measure') && isdineventmeasureerror && issynccondoff
        combinetiming('CondOff','Sync_SUFICI',ex.DisplayLatency);
        condoffversion='Sync';
    end
    if isfield(ex.CondTest,'CondOn') && strcmp(condonversion,'Sync') && isdineventsyncerror && iscommandcondon
        combinetiming('CondOn','Command_COND',ex.DisplayLatency);
    end
    if isfield(ex.CondTest,'CondOff') && strcmp(condoffversion,'Sync') && isdineventsyncerror && iscommandcondoff
        combinetiming('CondOff','Command_SUFICI',ex.DisplayLatency);
    end
    % Try to get first Cond On/Off timing
    if isfield(ex.CondTest,'CondOn')
        ex.CondTest.CondOn=firsteventtime('CondOn');
    end
    if isfield(ex.CondTest,'CondOff')
        ex.CondTest.CondOff=firsteventtime('CondOff');
    end
    % Parse CondOff when no PreICI/SufICI
    if isfield(ex.CondTest,'CondOn') && ~isfield(ex.CondTest,'CondOff')
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
    function ts = firsteventtime(event)
        ts=[];
        for i=1:nct
            vs=ex.CondTest.(event){i};
            if ~isempty(vs)
                t=vs(1);
            else
                t=NaN;
            end
            ts=[ts,t];
        end
    end
%% Try parse Environment Parameter
if ~isempty(ex.EnvParam)
    ex.EnvParam = tryparseparamstruct(ex.EnvParam);
end
%% Try parse Experiment Parameter
if ~isempty(ex.Param)
    ex.Param = tryparseparamstruct(ex.Param);
end
%% Parse CondTest Condition
isenvori=containsparam(ex.EnvParam,'Ori');
isenvorioffset = containsparam(ex.EnvParam,'OriOffset');
isenvposition=containsparam(ex.EnvParam,'Position');
isenvpositionoffset=containsparam(ex.EnvParam,'PositionOffset');
if containsparam(ex.EnvParam,'OriPositionOffset')
    isenvoripositionoffset=getparam(ex.EnvParam,'OriPositionOffset');
else
    isenvoripositionoffset=false;
end
if ~isempty(ex.Cond) && isfield(ex.CondTest,'CondIndex')
    % parse condition factor values
    fs = fieldnames(ex.Cond);
    for i=1:length(fs)
        f=fs{i};
        ex.Cond.(f) = cellfun(@(x)tryparseparam(f,x),ex.Cond.(f),'uniformoutput',false);
    end
    isori = ismember('Ori',fs);
    isorioffset = ismember('OriOffset',fs);
    isposition = ismember('Position',fs);
    ispositionoffset=ismember('PositionOffset',fs);
    % Initialize condtest condition
    for fi=1:length(fs)
        f=fs{fi};
        ctc.(f) = cell(1, nct);
    end
    if isori || isorioffset
        ctc.Ori_Final = cell(1, nct);
    end
    if isposition || ispositionoffset
        ctc.Position_Final = cell(1, nct);
    end
    % parse condition for each condtest
    for i=1:nct
        if ex.CondTest.CondIndex(i) < 1
            continue;
        end
        
        for fi=1:length(fs)
            f=fs{fi};
            ctc.(f)(i)=ex.Cond.(f)(ex.CondTest.CondIndex(i));
        end
        % parse final orientation
        if isori
            ori = ctc.Ori{i};
        elseif isenvori
            ori = getparam(ex.EnvParam,'Ori');
        else
            ori=0;
        end
        if isorioffset
            orioffset = ctc.OriOffset{i};
        elseif isenvorioffset
            orioffset = getparam(ex.EnvParam,'OriOffset');
        else
            orioffset=0;
        end
        if isori || isorioffset
            ctc.Ori_Final{i}=ori+orioffset;
        end
        % parse final position
        if isposition
            position=ctc.Position{i};
        elseif isenvposition
            position=getparam(ex.EnvParam,'Position');
        else
            position=zeros(1,3);
        end
        if ispositionoffset
            positionoffset=ctc.PositionOffset{i};
        elseif isenvpositionoffset
            positionoffset = getparam(ex.EnvParam,'PositionOffset');
        else
            positionoffset=zeros(1,3);
        end
        if isposition || ispositionoffset
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