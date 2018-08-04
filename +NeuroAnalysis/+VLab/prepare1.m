function [ex] = prepare1(ex,dataset)
%PREPARE1 Prepare version 1 of VLab data format
%   Detailed explanation goes here

import NeuroAnalysis.Base.* NeuroAnalysis.VLab.*
%% Trim CondTest
ctnames = fieldnames(ex.CondTest);
nct = min(cellfun(@(x)length(ex.CondTest.(x)),ctnames));
for i=1:length(ctnames)
    ex.CondTest.(ctnames{i}) = ex.CondTest.(ctnames{i})(1:nct);
end
%% Parse t0
ex.t0=0;
if ~isempty(dataset) && isfield(dataset,'digital')
    startsyncdchidx = find(arrayfun(@(x)x.channel==vlabconfig.StartSyncDCh,dataset.digital));
    if ~isempty(startsyncdchidx)
        ex.t0=dataset.digital(startsyncdchidx).time;
    end
end
%% Parse CondTest
if isfield(ex.CondTest,'CondIndex')
    % Convert to 1-base
    ex.CondTest.CondIndex =cellfun(@(x)int32(x+1), ex.CondTest.CondIndex);
end
if isfield(ex.CondTest,'CondRepeat')
    ex.CondTest.CondRepeat = cellfun(@(x)int32(x), ex.CondTest.CondRepeat);
end
if isfield(ex.CondTest,'BlockIndex')
    % Convert to 1-base
    ex.CondTest.BlockIndex =cellfun(@(x)int32(x+1), ex.CondTest.BlockIndex);
end
if isfield(ex.CondTest,'BlockRepeat')
    ex.CondTest.BlockRepeat = cellfun(@(x)int32(x), ex.CondTest.BlockRepeat);
end
    function searchrecover(from,to,data)
        names = fieldnames(ex.CondTest);
        fromnames = names(cellfun(@(x)startswith(x,from),names));
        for c=1:length(fromnames)
            e= replace(c,from,to);
            for i=1:nct
                recovered=[];
                ctses = ex.CondTest.(fromnames{c}){i};
                if ~isempty(ctses)
                    for j=1:length(ctses)
                        t=trysearchtime(ctses(j),data,vlabconfig.MaxDisplayLatencyError);
                        recovered=[recovered,t];
                    end
                end
                if ~isempty(recovered) && all(arrayfun(@(x)isnan(x),recovered))
                    recovered=[];
                end
                ex.CondTest.(e){i}=recovered;
            end
        end
    end
    function combinetiming(event,combineevent)
        for i=1:nct
            ctts=ex.CondTest.(event){i};
            cttcs = ex.CondTest.(combineevent){i};
            if isempty(ctts)
                if ~isempty(cttcs)
                    ex.CondTest.(event){i}=cttcs;
                end
            else
                if ~isempty(cttcs)
                    for j=1:length(ctts)
                        if isnan(ctts(j)) && ~isnan(cttcs(j))
                            ex.CondTest.(event){i}(j)=cttcs(j);
                        end
                    end
                end
            end
        end
    end
if isfield(ex.CondTest,'Event') && isfield(ex.CondTest,'SyncEvent')
    % Parse Sync Event VLab Timing
    ses=[];sectidx=[];isadddisplaylatency=true;
    for i = 1:nct
        ctses = ex.CondTest.SyncEvent{i};
        for j = 1:length(ctses)
            e = ctses{j};
            ses=[ses,{e}];
            sectidx=[sectidx,i];
            ex.CondTest.(['VLab_',e]){i} = arrayfun(@(x)vlabtime2reftime(ex,x,isadddisplaylatency),findeventtime(ex.CondTest.Event{i},e));
        end
    end
    % Check Sync Event Data
    if ~isempty(dataset) && isfield(dataset,'digital')
        if ex.EventSyncProtocol.nSyncChannel==1 && ex.EventSyncProtocol.nSyncpEvent==1
            eventsyncdchidx = find(arrayfun(@(x)x.channel==vlabconfig.EventSyncDCh,dataset.digital));
            eventmeasuredchidx = find(arrayfun(@(x)x.channel==vlabconfig.EventMeasureDCh,dataset.digital));
            
            isdineventsync = ~isempty(eventsyncdchidx);
            if isdineventsync
                dineventsynctime = dataset.digital(eventsyncdchidx).time;
                dineventsyncdata = dataset.digital(eventsyncdchidx).data;
                if length(ses)==length(dineventsyncdata) && all(diff(dineventsyncdata))
                    isdineventsyncerror=false;
                else
                    isdineventsyncerror=true;
                end
            end
            
            isdineventmeasure = ~isempty(eventmeasuredchidx);
            if isdineventmeasure
                dineventmeasuretime = dataset.digital(eventmeasuredchidx).time;
                dineventmeasuredata = dataset.digital(eventmeasuredchidx).data;
                if length(ses)==length(dineventmeasuredata) && all(diff(dineventmeasuredata))
                    isdineventmeasureerror=false;
                else
                    isdineventmeasureerror=true;
                end
            end
        end
    end
    % Parse Sync Event Timing
    if isdineventsync
        if ~isdineventsyncerror
            if isadddisplaylatency
                displaylatency=ex.DisplayLatency;
            else
                displaylatency=0;
            end
            for i=1:length(ses)
                se = ['Sync_',ses{i}];
                if ~isfield(ex.CondTest,se)
                    ex.CondTest.(se)=cell(1,nct);
                end
                ex.CondTest.(se){sectidx(i)} = [ex.CondTest.(se){sectidx(i)},dineventsynctime(i)+displaylatency];
            end
        else
            % Try to recover as many sync timing as possible based on VLab Timing
            searchrecover('VLab_','Sync_',dineventsynctime);
        end
    end
    % Parse Sync Event Measure Timing
    if isdineventmeasure
        if ~isdineventmeasureerror
            for i=1:length(ses)
                se = ['Measure_',ses{i}];
                if ~isfield(ex.CondTest,se)
                    ex.CondTest.(se)=cell(1,nct);
                end
                ex.CondTest.(se){sectidx(i)} = [ex.CondTest.(se){sectidx(i)},dineventmeasuretime(i)];
            end
        else
            % Try to recover as many measure timing as possible based on Sync Timing
            searchrecover('Sync_','Measure_',dineventmeasuretime);
        end
    end
    % Try to get the most accurate and complete Cond On/Off Time
    ismeasurecondon = isfield(ex.CondTest,'Measure_COND');
    issynccondon = isfield(ex.CondTest,'Sync_COND');
    isvlabcondon = isfield(ex.CondTest,'VLab_COND');
    ismeasurecondoff = isfield(ex.CondTest,'Measure_SUFICI');
    issynccondoff = isfield(ex.CondTest,'Sync_SUFICI');
    isvlabcondoff = isfield(ex.CondTest,'VLab_SUFICI');
    
    condonversion='None';condoffversion='None';
    if ismeasurecondon
        ex.CondTest.CondOn=ex.CondTest.Measure_COND;
        condonversion='Measure';
    elseif issynccondon
        ex.CondTest.CondOn=ex.CondTest.Sync_COND;
        condonversion='Sync';
    elseif isvlabcondon
        ex.CondTest.CondOn=ex.CondTest.VLab_COND;
    end
    if ismeasurecondoff
        ex.CondTest.CondOff=ex.CondTest.Measure_SUFICI;
        condoffversion='Measure';
    elseif issynccondoff
        ex.CondTest.CondOff=ex.CondTest.Sync_SUFICI;
        condoffversion='Sync';
    elseif isvlabcondoff
        ex.CondTest.CondOff=ex.CondTest.VLab_SUFICI;
    end
    
    if isfield(ex.CondTest,'CondOn') && strcmp(condonversion,'Measure') && isdineventmeasureerror && issynccondon
        combinetiming('CondOn','Sync_COND');
        condonversion='Sync';
    end
    if isfield(ex.CondTest,'CondOff') && strcmp(condoffversion,'Measure') && isdineventmeasureerror && issynccondoff
        combinetiming('CondOff','Sync_SUFICI');
        condoffversion='Sync';
    end
    if isfield(ex.CondTest,'CondOn') && strcmp(condonversion,'Sync') && isdineventsyncerror && isvlabcondon
        combinetiming('CondOn','VLab_COND');
    end
    if isfield(ex.CondTest,'CondOff') && strcmp(condoffversion,'Sync') && isdineventsyncerror && isvlabcondoff
        combinetiming('CondOff','VLab_SUFICI');
    end
    % Parse CondOff when no PreICI/SufICI
    if isfield(ex.CondTest,'CondOn') && ~isfield(ex.CondTest,'CondOff')
        for i=1:nct-1
            currentontime=ex.CondTest.CondOn{i};
            nextontime = ex.CondTest.CondOn{i+1};
            if (nextontime - currentontime) > (ex.CondDur+2*vlabconfig.MaxDisplayLatencyError)
                ex.CondTest.CondOff{i} = currentontime + ex.CondDur;
            else
                ex.CondTest.CondOff{i}=nextontime;
            end
        end
        ex.CondTest.CondOff{nct}=ex.CondTest.CondOn{nct}+ex.CondDur;
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
    % parse condition for each condtest
    for i=1:nct
        for fi=1:length(fs)
            f=fs{fi};
            ctc.(f)(i)=ex.Cond.(f)(ex.CondTest.CondIndex(i));
        end
        % parse final orientation
        isori = ismember('Ori',fs);
        isorioffset = ismember('OriOffset',fs);
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
        isposition = ismember('Position',fs);
        ispositionoffset=ismember('PositionOffset',fs);
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