function [ex] = parseex(ex)
%PARSEEX Parse VLab Experimental Data
%   Detailed explanation goes here

import NeuroAnalysis.VLab.*
% Trim CondTest
ctnames = fieldnames(ex.CondTest);
nct = min(cellfun(@(x)length(ex.CondTest.(x)),ctnames));
for i=1:length(ctnames)
    ex.CondTest.(ctnames{i}) = ex.CondTest.(ctnames{i})(1:nct);
end
% Convert CondIndex to 1-base
ex.CondTest.CondIndex =cellfun(@(x)int32(x+1), ex.CondTest.CondIndex);
ex.CondTest.CondRepeat = cellfun(@(x)int32(x), ex.CondTest.CondRepeat);
% Try parse Environment Parameter
if ~isempty(ex.EnvParam)
    envparamnames = fieldnames(ex.EnvParam);
    for i=1:length(envparamnames)
        p = envparamnames{i};
        ex.EnvParam.(p) = tryparseparam(p,ex.EnvParam.(p));
    end
end
% Try parse User Experiment Parameter
if ~isempty(ex.Param)
    paramnames = fieldnames(ex.Param);
    for i=1:length(paramnames)
        p = paramnames{i};
        ex.Param.(p) = tryparseparam(p,ex.Param.(p));
    end
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
    % try parse condition factor value
    fs = fieldnames(ex.Cond);
    for i=1:length(fs)
        f=fs{i};
        ex.Cond.(f) = cellfun(@(x)tryparseparam(f,x),ex.Cond.(f),'uniformoutput',0);
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

ex.CondTest.PreICIOnTime =cellfun(@(x)vlabtimetodatatime(ex,findstatetime(x,'PREICI'),true),ex.CondTest.CONDSTATE);
ex.CondTest.CondOnTime = cellfun(@(x)vlabtimetodatatime(ex,findstatetime(x,'COND'),true),ex.CondTest.CONDSTATE);
ex.CondTest.SufICIOnTime = cellfun(@(x)vlabtimetodatatime(ex,findstatetime(x,'SUFICI'),true),ex.CondTest.CONDSTATE);
end

