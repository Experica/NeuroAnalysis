function [ex] = parseex(ex)
%PARSEEX Parse version 0 of Experica.Command Data
%   Detailed explanation goes here

import NeuroAnalysis.Base.* NeuroAnalysis.Experica.*
% Trim CondTest
ctnames = fieldnames(ex.CondTest);
nct = min(cellfun(@(x)length(ex.CondTest.(x)),ctnames));
for i=1:length(ctnames)
    ex.CondTest.(ctnames{i}) = ex.CondTest.(ctnames{i})(1:nct);
end

    function [t] = findstatetime(stateevents,state)
        %FINDSTATETIME Extract timestamps of state event
        
        t=[];
        for i=1:length(stateevents)
            if isfield(stateevents{i},state)
                t = [t,stateevents{i}.(state)];
            end
        end
        if isempty(t)
            t=NaN;
        end
        
    end

    function [datatime] = todatatime(ex,time,isaddlatency)
        %TODATATIME Convert to Data Timing
        
        if nargin==2
            isaddlatency=false;
        end
        
        datatime = time*(1+ex.TimerDriftSpeed)+ex.t0;
        if isaddlatency
            datatime = datatime + ex.Latency;
        end
        
    end

ex.CondTest.CondIndex =cellfun(@(x)int32(x+1), ex.CondTest.CondIndex); % Convert to 1-base
ex.CondTest.CondRepeat = cellfun(@(x)int32(x), ex.CondTest.CondRepeat);
ex.CondTest.PreICIOnTime =cellfun(@(x)todatatime(ex,findstatetime(x,'PREICI'),true),ex.CondTest.CONDSTATE);
ex.CondTest.CondOnTime = cellfun(@(x)todatatime(ex,findstatetime(x,'COND'),true),ex.CondTest.CONDSTATE);
ex.CondTest.SufICIOnTime = cellfun(@(x)todatatime(ex,findstatetime(x,'SUFICI'),true),ex.CondTest.CONDSTATE);
if isfield(ex.CondTest,'TRIALSTATE')
    ex.CondTest.PreITIOnTime =cellfun(@(x)todatatime(ex,findstatetime(x,'PREITI'),true),ex.CondTest.TRIALSTATE);
    ex.CondTest.TrialOnTime = cellfun(@(x)todatatime(ex,findstatetime(x,'TRIAL'),true),ex.CondTest.TRIALSTATE);
    ex.CondTest.SufITIOnTime = cellfun(@(x)todatatime(ex,findstatetime(x,'SUFITI'),true),ex.CondTest.TRIALSTATE);
end
if isfield(ex.CondTest,'BLOCKSTATE')
    ex.CondTest.PreIBIOnTime =cellfun(@(x)todatatime(ex,findstatetime(x,'PREIBI'),true),ex.CondTest.BLOCKSTATE);
    ex.CondTest.BlockOnTime = cellfun(@(x)todatatime(ex,findstatetime(x,'BLOCK'),true),ex.CondTest.BLOCKSTATE);
    ex.CondTest.SufIBIOnTime = cellfun(@(x)todatatime(ex,findstatetime(x,'SUFIBI'),true),ex.CondTest.BLOCKSTATE);
end

% Try parse Environment Parameter
if ~isempty(ex.EnvParam)
    ex.EnvParam = tryparseparamstruct(ex.EnvParam);
end
% Try parse Experiment Parameter
if ~isempty(ex.Param)
    ex.Param = tryparseparamstruct(ex.Param);
end

isenvori=containsparam(ex.EnvParam,'Ori');
isenvorioffset = containsparam(ex.EnvParam,'OriOffset');
isenvposition=containsparam(ex.EnvParam,'Position');
isenvpositionoffset=containsparam(ex.EnvParam,'PositionOffset');
if containsparam(ex.EnvParam,'OriPositionOffset')
    isenvoripositionoffset=getparam(ex.EnvParam,'OriPositionOffset');
else
    isenvoripositionoffset=false;
end
if ~isempty(ex.Cond)
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

end