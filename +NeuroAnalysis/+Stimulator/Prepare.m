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
%% check analyzer and log file
stimulatordataset = struct([]);
islogfile=false;
[filedir,filename,ext] = fileparts(filepath);
fns = strsplit(filename,'_');
fnu = fns{2};
fns{2}=fnu(2:end);
logfilepath = fullfile(filedir,['*_',strjoin(fns,'_'),'.mat']);
d = dir(logfilepath);
if ~isempty(d)
    islogfile=true;
    logfilepath = fullfile(d(1).folder,d(1).name);
end

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
if islogfile
    disp(['Reading Stimulator Log File:    ',logfilepath,'    ...']);
    ex.raw.log = load(logfilepath,'-mat');
    disp(['Reading Stimulator Log File:    ',logfilepath,'    Done.']);
end
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
        %% Parse Stimulator Analyzer
        
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
        ex.EnvParam.ScreenSizeDegree = [deg*ex.EnvParam.ScreenAspect*2, deg*2];
        ex.EnvParam.Size = [ex.EnvParam.x_size,ex.EnvParam.y_size];
        
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
        ex.nTrial = length(cti);
        ex.CondTestCond=ctc;
        %% Parse Stimulator Log
        if isfield(ex.raw,'log')
            log = ex.raw.log;
            ex.EnvParam.DisplayRefreshRate = log.frate;
            if strcmp(ex.ID,'FG')
                ts = fieldnames(log);
                ts=ts(cellfun(@(x)startsWith(x,'randlog'),ts));
                tn = cellfun(@(x)str2double(x(10:end)),ts);
                [~,ti] = sort(tn);
                ts=ts(ti);
                
                condidx=[];
                ctc=[];
                cond=[];
                for t=1:length(ts)
                    tp = log.(ts{t});
                    if isempty(cond)
                        for i =1:size(tp.domains.Cond,1)
                            cond=[cond,parsehartley(tp.domains.Cond(i,:),ex.EnvParam.Size)];
                        end
                    end
                    for c=1:length(tp.seqs.frameseq)
                        ci = tp.seqs.frameseq(c);
                        condidx=[condidx,ci];
                        ctc=[ctc,cond(ci)];
                    end
                end
                ex.CondTest.CondIndex=int64(condidx);
                ex.CondTestCond = NeuroAnalysis.Base.arraystruct2structarray(ctc);
            end
        end
        
        %% Convert Grating Drifting Direction to Orientation
        if strcmp(ex.ID,'PG')
            if isfield(ex.CondTestCond,'Ori')
                ex.CondTestCond.Ori = mod(ex.CondTestCond.Ori+270,360);
            end
        end
        %%
        function [c]=parsehartley(p,size)
            oridom = p(1);kx=p(2);ky=p(3);bwdom=p(4);colordom=p(5);
            if kx==0 && ky~=0   % 0
                c.Ori=0;
                sf = abs(ky)/size(2);
            elseif ky==0 && kx~=0   % 90
                c.Ori=90;
                sf = abs(kx)/size(1);
            elseif ky*ky > 0   % (90 180)
                c.Ori = 180 - atand(ky/kx);
                sf = abs(ky/(size(2)*sin(c.Ori-90)))
            elseif kx*ky < 0   % (0 90)
                c.Ori = abs(atand(ky/kx));
                sf = abs(kx/(size(1)*sin(c.Ori)))
            end

            % if 0 < c.Ori <= 45 || 135 <= c.Ori < 180
            %     sf = sqrt(kx^2 + ky^2) / abs((size(1)/cos(c.Ori));
            % elseif 45 < c.Ori < 90 || 90 < c.Ori < 135
            %     sf = sqrt(kx^2 + ky^2) / (size(2)/sin(c.Ori))
            c.SpatialFreq = sf;

            if ky>=0
                q=0;
                if kx<0 && ky==0
                    q=0.25;
                end
            else
                q=0.25;
            end
            if bwdom==-1
                q=q+0.5;
            end
            c.SpatialPhase = q;
        end
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
                
                if ex.nTrial == length(ex.nidq.digital(1).time)/2
                    if strcmp(ex.ID,'FG') && ex.PreICI==0 && ex.SufICI==0
                        ex.CondTest.TrialOn = sample2time(ex.nidq.digital(1).time(1:2:end),ex.nidq.fs,dataset.secondperunit)+ex.PreITI;
                        ex.CondTest.TrialOff= sample2time(ex.nidq.digital(1).time(2:2:end),ex.nidq.fs,dataset.secondperunit)-ex.SufITI;
                        odt=ex.nidq.digital(2).time;
                        
                        % the first oneframepulse of trial will merge with the first condition flip if PreITI==0, so there
                        % would be 2 flip instead of 4 flip of trial oneframepulse
                        if ex.PreITI==0
                            ncondintrial = (length(odt)-2*ex.nTrial)/ex.nTrial;
                            condon=[];
                            for i=1:ex.nTrial
                                t = odt((1:ncondintrial)+(i-1)*(ncondintrial+2));
                                t(1)=t(1)+pd;
                                condon=[condon,t];
                            end
                            ex.CondTest.CondOn = sample2time(condon,ex.nidq.fs,dataset.secondperunit);
                        else
                            oneframepulseindex=find(diff(odt)<pd*1.3);
                            odt([oneframepulseindex,oneframepulseindex+1])=[];
                            ex.CondTest.CondOn = sample2time(odt,ex.nidq.fs,dataset.secondperunit);
                        end
                        
                        if isfield(ex.EnvParam,'DisplayRefreshRate')
                            ex.CondDur = ex.CondDur/ex.EnvParam.DisplayRefreshRate;
                        else
                            ex.CondDur = ex.CondDur*pd/ex.nidq.fs;
                        end
                        ex.CondTest.CondOff= [ex.CondTest.CondOn(2:end),ex.CondTest.CondOn(end)+ex.CondDur];
                        iidx = (ex.CondTest.CondOff-ex.CondTest.CondOn) > pd/ex.nidq.fs+ex.CondDur;
                        if any(iidx)
                            ex.CondTest.CondOff(iidx) = ex.CondTest.CondOn(iidx)+ex.CondDur;
                        end
                        iidx = (ex.CondTest.CondOff-ex.CondTest.CondOn) < ex.CondDur-pd/ex.nidq.fs;
                        if any(iidx)
                            ex.CondTest.CondOff(iidx) = ex.CondTest.CondOn(iidx)+ex.CondDur;
                        end
                    else
                        if ex.nTrial == length(ex.nidq.digital(2).time)/4 % photodiode
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