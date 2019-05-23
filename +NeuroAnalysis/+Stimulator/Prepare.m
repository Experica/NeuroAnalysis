function [ stimulatordataset ] = Prepare( filepath,varargin )
%PREPARE Read and Prepare Stimulator experimental data
%   Detailed explanation goes here

p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addOptional(p,'dataset',struct([]));
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
dataset = p.Results.dataset;
%% check file
stimulatordataset = struct([]);
if ~exist(filepath, 'file')
    error(['Analyzer file: ',filepath,' not found.']);
end
d = dir(filepath);
%% Read data
disp(['Reading Stimulator File:    ',filepath,'    ...']);
ex = struct;
ex.raw = load(filepath,'-mat');
ex.raw = ex.raw.Analyzer;
ex.source = filepath;
ex.sourceformat = 'Stimulator';
disp(['Reading Stimulator File:    ',filepath,'    Done.']);
%% Prepare data
disp(['Preparing Stimulator File:    ',filepath,'    ...']);
if isempty(ex.raw)
    warning(['Preparing Stimulator File:    ',filepath,'    Empty.']);
    return;
end
ex = parseex( ex );
if ~isempty(dataset) && (isfield(dataset,'ap') || isfield(dataset,'lf') || isfield(dataset,'nidq'))
    ex = parseevent(ex);
end
% Standardize experimental parameters
ex = NeuroAnalysis.Base.StandardizeEx(ex);
% Organize into dataset
stimulatordataset = ex;
stimulatordataset.date = d.datenum;

disp(['Preparing Stimulator File:    ',filepath,'    Done.']);
%%
    function [ ex ] = parseex( ex )
        % Parse Stimulator Analyzer
        
        m = ex.raw.M;
        p = ex.raw.P;
        l = ex.raw.L;
        lp= ex.raw.loops;
        % Experimental parameters
        ex.Subject_ID = m.anim;
        ex.RecordSession = '';
        ex.RecordSite = m.unit;
        ex.Hemisphere = m.hemi;
        ex.DisplayType= m.monitor;
        ex.ID = p.type;
        
        if l.rand
            ex.CondSampling = 'UniformWithoutReplacement';
        else
            ex.CondSampling = 'Ascending';
        end
        ex.BlockSampling = 'Ascending';
        ex.CondRepeat = l.reps;
        ex.BlockRepeat = 1;
        ex.BlockParam = [];
        
        for i=1:length(l.param)
            ip = l.param{i};
            ex.Factor.(ip{1})=ip{2};
        end
        
        % Environmental parameters
        for i=1:length(p.param)
            ip=p.param{i};
            ex.EnvParam.(ip{1})=ip{3};
        end
        switch ex.ID
            case 'FG'
                ex.PreITI = ex.EnvParam.predelay;
                ex.TrialDur = ex.EnvParam.stim_time;
                ex.SufITI = ex.EnvParam.postdelay;
                ex.PreICI = 0;
                ex.CondDur = ex.EnvParam.h_per;
                ex.SufICI = 0;
            case 'PG'
                ex.PreICI = ex.EnvParam.predelay;
                ex.CondDur = ex.EnvParam.stim_time;
                ex.SufICI = ex.EnvParam.postdelay;
        end
        
        ex.EnvParam.ScreenToEye = m.screenDist;
        ex.EnvParam.ScreenResolution = [m.xpixels m.ypixels];
        ex.EnvParam.ScreenSize = [m.screenXcm m.screenYcm];
        ex.EnvParam.MarkerSize = m.syncSize;
        
        ex.EnvParam.ScreenAspect = ex.EnvParam.ScreenResolution(1)/ex.EnvParam.ScreenResolution(2);
        ex.EnvParam.ScreenHalfHeight = m.screenYcm/2;
        deg = atand(ex.EnvParam.ScreenHalfHeight/ex.EnvParam.ScreenToEye);
        ex.EnvParam.ScreenDegrees = [deg*ex.EnvParam.ScreenAspect*2, deg*2];
        
        % Condition Tests
        ctc = struct;
        cti=[];
        factornames = fieldnames(ex.Factor);
        for i=1:length(factornames)
            ctc.(factornames{i})={};
        end
        for i=1:length(lp.conds)
            ic = lp.conds{i};
            for j=1:length(ic.repeats)
                cti=[cti,ic.repeats{j}.trialno];
                for k=1:length(ic.symbol)
                    sym = ic.symbol{k};
                    val = ic.val{k};
                    if strcmp(sym,'blank')
                        sym=factornames{k};
                        val='blank';
                    end
                    ctc.(sym)=[ctc.(sym),val];
                end
            end
        end
        [~,si]=sort(cti);
        for i=1:length(factornames)
            ctc.(factornames{i})=ctc.(factornames{i})(si);
        end
        ex.nCondTest = length(cti);
        ex.CondTestCond=ctc;
        
    end
%%
    function [ex] = parseevent(ex)
        % Parse Experimental Event Time synced and measured in each data stream
        import NeuroAnalysis.Base.*
        
        if isfield(dataset,'nidq') && ~isempty(dataset.nidq.meta.fileName)
            nsample=double(dataset.nidq.meta.nFileSamp);
            fs = dataset.nidq.meta.fs;
            chn = double(dataset.nidq.meta.nSavedChans);
            if chn>0
                binmap = memmapfile(dataset.nidq.meta.fileName,'Format',{'int16',[chn,nsample],'nidq'});
                
                ex.nidq.fs=fs;
                [di,dv]=parsedigitalinanalog(binmap.Data.nidq(2,:),nsample,2520,1080);
                ex.nidq.digital(1).channel=1;
                ex.nidq.digital(1).time=di;
                ex.nidq.digital(1).value=dv;
                
                if strcmpi(ex.DisplayType,'crt')
                    [di,dv,pd]=parsecrtdigitalinanalog(binmap.Data.nidq(1,:),300,1500,1300);
                    % correct CRT scan, because upperleft corner photodiode only detects the first scan line.
                    di = di + round(pd*ex.EnvParam.y_pos/ex.EnvParam.ScreenResolution(2));
                else
                    [di,dv]=parsedigitalinanalog(binmap.Data.nidq(1,:),nsample,1500,1300);
                end
                ex.nidq.digital(2).channel=2;
                ex.nidq.digital(2).time=di;
                ex.nidq.digital(2).value=dv;
                
                if ex.nCondTest == length(ex.nidq.digital(1).time)/2
                    if strcmp(ex.ID,'FG') && ex.PreICI==0 && ex.SufICI==0
                        ex.CondTest.TrialOn = sample2time(ex.nidq.digital(1).time(1:2:end),ex.nidq.fs,dataset.secondperunit)+ex.PreITI;
                        ex.CondTest.TrialOff= sample2time(ex.nidq.digital(1).time(2:2:end),ex.nidq.fs,dataset.secondperunit)-ex.SufITI;
                        odt=ex.nidq.digital(2).time;
                        oneframepulseindex=find(diff(odt)<pd*1.5);
                        odt([oneframepulseindex,oneframepulseindex+1])=[];
                        ex.CondTest.CondOn = sample2time(odt,ex.nidq.fs,dataset.secondperunit);
                        ex.CondDur = ex.CondDur*pd/ex.nidq.fs;
                        ex.CondTest.CondOff= [ex.CondTest.CondOn(2:end),ex.CondTest.CondOn(end)+ex.CondDur];
                    else
                        if ex.nCondTest == length(ex.nidq.digital(2).time)/4 % photodiode
                            ex.CondTest.CondOn=sample2time(ex.nidq.digital(2).time(1:4:end),ex.nidq.fs,dataset.secondperunit)+ex.PreICI;
                            ex.CondTest.CondOff=sample2time(ex.nidq.digital(2).time(3:4:end),ex.nidq.fs,dataset.secondperunit);
                        else % digital
                            ex.CondTest.CondOn = sample2time(ex.nidq.digital(1).time(1:2:end),ex.nidq.fs,dataset.secondperunit)+ex.PreICI;
                            ex.CondTest.CondOff= sample2time(ex.nidq.digital(1).time(2:2:end),ex.nidq.fs,dataset.secondperunit)-ex.SufICI;
                        end
                    end
                end
            end
        end
        
        if isfield(dataset,'lf') && ~isempty(dataset.lf.meta.fileName)
            nsample=double(dataset.lf.meta.nFileSamp);
            fs = dataset.lf.meta.fs;
            chn = double(dataset.lf.meta.nSavedChans);
            if dataset.lf.meta.snsApLfSy(3)>0
                binmap = memmapfile(dataset.lf.meta.fileName,'Format',{'uint16',[chn,nsample],'lf'});
                ex.lf.digital=parsedigitalbitinanalog(binmap.Data.lf(chn,:),nsample,16);
                ex.lf.fs=fs;
            end
        end
        
        % if isfield(dataset,'ap') && ~isempty(dataset.ap.meta.fileName)
        %     nsample=dataset.ap.meta.nFileSamp;
        %     fs = dataset.ap.meta.imSampRate;
        %     chns = dataset.ap.meta.snsApLfSy;
        %     chn = dataset.ap.meta.nSavedChans;
        %     if chns(3)>0
        %         binmap = memmapfile(dataset.ap.meta.fileName,'Format',{'uint16',[chn,nsample],'ap'});
        %         ex.ap.digital=parsedigital(binmap.Data.ap(chn,:),nsample,16);
        %         ex.ap.fs=fs;
        %     end
        % end
        
    end
end