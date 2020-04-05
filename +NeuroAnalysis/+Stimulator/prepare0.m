function [ stimulatordataset ] = prepare0( filepath,varargin )
% function [ ex ] = prepare0( ex,dataset )
%PREPARE Read and Prepare Stimulator experimental data
%   Detailed explanation goes here
% !!!Only works with old 2P data (Before AF5)!!! Peichao

p = inputParser;
p.StructExpand = false;
addRequired(p,'filepath');
addOptional(p,'dataset',struct([]));
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
dataset = p.Results.dataset{1,1};
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

% Parse spike data
if ~isempty(dataset) && strcmp(dataset.exptId,'Hartley') && strcmp(dataset.eventSource,'spike2')
    ex = parseevent(ex);
end

% Align scanbox frame with spike2 for Hartley only 
if ~isempty(dataset) && strcmp(dataset.exptId,'Hartley') && strcmp(dataset.eventSource,'spike2')
    ex = parsescanbox(ex);
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
        ex.RecordSite = ['u',m.unit];
        ex.Hemisphere = m.hemi;
        ex.TestID = m.expt;
        ex.DisplayType= m.monitor;
        ex.ID = p.type;
        ex.Name = ex.ID;
        
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
            if strcmp(ex.ID, 'PG') && strcmp(ip{1}, 'ori')
                ip{1} = 'Dir';
            end
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
                ex.CondDur = ex.EnvParam.h_per;   % in frame number, not in time, PL
                ex.SufICI = 0;
            case 'PG'
                ex.PreICI = ex.EnvParam.predelay;
                ex.CondDur = ex.EnvParam.stim_time;  % here is in time, sec
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
        for i=1:length(lp.conds)   % All conditions test in expt
            ic = lp.conds{i};
            for j=1:length(ic.repeats)
                cti=[cti,ic.repeats{j}.trialno];
                for k=1:length(ic.symbol)
                    sym = ic.symbol{k};
                    if strcmp(ex.ID, 'PG') && strcmp(sym, 'ori')
                        sym = 'Dir';
                    end
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
        if isfield(ex.CondTestCond,'colormod') && strcmp(ex.ID, 'PG')
            ex.ID = 'DirSFColor';
        elseif ~isfield(ex.CondTestCond,'colormod') && strcmp(ex.ID, 'PG')
            ex.ID = 'DirSF';
        end
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
                ex.Cond = NeuroAnalysis.Base.arraystruct2structarray(cond);
                ex.CondTest.CondIndex=int64(condidx);
                ex.CondTestCond = NeuroAnalysis.Base.arraystruct2structarray(ctc);
            end
        end
        
        %% Convert Grating Drifting Direction to Orientation
%         if strcmp(ex.ID,'PG')
%             if isfield(ex.CondTestCond,'Dir')
%                 ex.CondTestCond.Ori = mod(ex.CondTestCond.Dir-90,180);
%             end
%         end
        %% Get Ori, SpatialFreq and SpatialPhase from Hartley Space parameters(horizontal ori = 0), PL
        function [c]=parsehartley(p,size)
            oridom = p(1);kx=p(2);ky=p(3);bwdom=p(4);colordom=p(5);
            if kx==0 && ky~=0   % 0
                c.Ori=0;
            elseif ky==0 && kx~=0   % 90
                c.Ori=90;
            elseif kx*ky > 0   % (0 90)
                c.Ori = atand(kx/ky); 
            elseif kx*ky < 0   % (90 180)
                c.Ori = 180 - atand(abs(kx/ky));
            else
                c.Ori=NaN;
            end
            c.SpatialFreq = sqrt((kx/size(1))^2 + (ky/size(2))^2);
            
            % The phase here is the phase calculated at the edge of image, PL
            % if kx>=0
            %     q=0.125;
            %     if ky<0 && kx==0
            %         q=0.375;
            %     end
            % else
            %     q=0.375;
            % end
            % if bwdom==-1
            %     q=q+0.5;
            % end

            % The phase here is the phase calculated at the center of image for regeneration in Julia & Experica, PL
            if mod(kx+kx,2) == 0
                if kx>=0
                    q=0.125;
                    if ky<0 && kx==0
                        q=0.375;
                    end
                else
                    q=0.375;
                end
                if bwdom==-1
                    q=q+0.5;
                end
            elseif mod(kx+kx,2) ~= 0
                if kx>=0
                    q=-0.375;
                    if ky<0 && kx==0
                        q=-0.125;
                    end
                else
                    q=-0.125;
                end
                if bwdom==-1
                    q=q-0.5;
                end
            end
            q=0.5-q; % phase 0 starts at pi(0.5) in Julia and Experica
            q = rem(q+1,1);  % transform negative to positive

            c.SpatialPhase = q;
        end
    end
%%
    function [ex] = parseevent(ex)
        % Parse Experimental Event Time synced and measured in each data stream
        import NeuroAnalysis.Base.*
        if isfield(dataset,'nidq')
            daqType = 'nidq';
        elseif strcmp(dataset.eventSource,'spike2')
            daqType = 'spike2';
        end
            
            
        if (isfield(dataset,'nidq') && ~isempty(dataset.nidq.meta.fileName)) || strcmp(dataset.eventSource,'spike2')
            switch daqType
                case 'nidq'
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
                    end
                case 'spike2'
                        nsample=size(dataset.spike2Data.time,1);
                        fs = dataset.spike2Data.fs;
                        chn = size(dataset.spike2Data.chanName,2);
                        digData = dataset.spike2Data.fs1208';  % Digital TTL
                        anaData = dataset.spike2Data.photodio'; % Photodiode
                        
                        if chn>0
                            ex.spike2.fs=fs;
                            [di,dv]=parsedigitalinanalog(digData,nsample,0.8,1.5);
                            ex.spike2.digital(1).channel=1;  % Digital
                            ex.spike2.digital(1).time=di;
                            ex.spike2.digital(1).value=dv;

                            if strcmpi(ex.DisplayType,'crt')
                                [di,dv,pd]=parsecrtdigitalinanalog(anaData,0.05,0.25,0.15);
                                % correct CRT scan, because upperleft corner photodiode only detects the first scan line.
                                di = di + round(pd*ex.EnvParam.y_pos/ex.EnvParam.ScreenResolution(2));
                            else
                                [di,dv]=parsedigitalinanalog(anaData,nsample,0.25,0.15);
                            end
                            ex.spike2.digital(2).channel=2;    % Photodiode
                            ex.spike2.digital(2).time=di;
                            ex.spike2.digital(2).value=dv;
%                             ex.spike2.digital(2).time=[di(1),di(1)+215, di(2:end)];
%                             ex.spike2.digital(2).value=[dv(1),dv(1),dv(2:end)];

                        end
            end
                
                if ex.nTrial == length(ex.(daqType).digital(1).time)/2
                    if strcmp(ex.ID,'FG') && ex.PreICI==0 && ex.SufICI==0
                        ex.CondTest.TrialOn = sample2time(ex.(daqType).digital(1).time(1:2:end),ex.(daqType).fs,dataset.secondperunit); %+ex.PreITI;
                        ex.CondTest.TrialOff= sample2time(ex.(daqType).digital(1).time(2:2:end),ex.(daqType).fs,dataset.secondperunit); %-ex.SufITI;
                        odt=ex.(daqType).digital(2).time;  % Photodiode channel time
                        
                        % the first oneframepulse of trial will merge with the first condition flip if PreITI==0, so there
                        % would be 2 flip instead of 4 flip of trial oneframepulse
                        if isfield(ex.raw,'log')
                            ex.ID = 'Hartley';
%                             ncondintrial =
%                             (length(odt)-2*ex.nTrial)/ex.nTrial;  % Peichao: This is for old Ephys only (before AF4)
%                             condon=[];
%                             for i=1:ex.nTrial
%                                 t = odt((1:ncondintrial)+(i-1)*(ncondintrial+2));
%                                 t(1)=t(1)+pd;
%                                 condon=[condon,t];
%                             end
%                             ex.CondTest.CondOn = sample2time(condon,ex.(daqType).fs,dataset.secondperunit);
                            % Peichao
                             nFlipPerTrial = length(odt)/ex.nTrial;  
                             odtTemp = reshape(odt, nFlipPerTrial, ex.nTrial);
                             odtTemp = odtTemp(3:end-2, :);  % Peichao: In old 2P data, photodiode channel has two more filps mark the start of pre and end of post.
                             condon= reshape(odtTemp, 1, []);
                             ex.CondTest.CondOn = sample2time(condon,ex.(daqType).fs,dataset.secondperunit);

                        else
                            ex.ID = 'Flash';
                            oneframepulseindex=find(diff(odt)<pd*1.3);
                            odt([oneframepulseindex,oneframepulseindex+1])=[];
                            ex.CondTest.CondOn = sample2time(odt,ex.(daqType).fs,dataset.secondperunit);
                        end
                        
                        if isfield(ex.EnvParam,'DisplayRefreshRate')
                            ex.CondDur = ex.CondDur/ex.EnvParam.DisplayRefreshRate;
                        else
                            ex.CondDur = ex.CondDur*pd/ex.(daqType).fs;
                        end
                        ex.CondTest.CondOff= [ex.CondTest.CondOn(2:end),ex.CondTest.CondOn(end)+ex.CondDur];
                        iidx = (ex.CondTest.CondOff-ex.CondTest.CondOn) > pd/ex.(daqType).fs+ex.CondDur;
                        if any(iidx)
                            ex.CondTest.CondOff(iidx) = ex.CondTest.CondOn(iidx)+ex.CondDur;
                        end
                        iidx = (ex.CondTest.CondOff-ex.CondTest.CondOn) < ex.CondDur-pd/ex.(daqType).fs;
                        if any(iidx)
                            ex.CondTest.CondOff(iidx) = ex.CondTest.CondOn(iidx)+ex.CondDur;
                        end
                    else
                        if isfield(ex.CondTestCond,'colormod')
                            ex.ID = 'DirSFColor';
                        else
                            ex.ID = 'DirSF';
                        end
                        if ex.nTrial == length(ex.(daqType).digital(2).time)/4 % photodiode
                            ex.CondTest.CondOn=sample2time(ex.(daqType).digital(2).time(1:4:end),ex.(daqType).fs,dataset.secondperunit)+ex.PreICI;
                            ex.CondTest.CondOff=sample2time(ex.(daqType).digital(2).time(3:4:end),ex.(daqType).fs,dataset.secondperunit);
                        else % digital
                            ex.CondTest.CondOn = sample2time(ex.(daqType).digital(1).time(1:2:end),ex.(daqType).fs,dataset.secondperunit)+ex.PreICI;
                            ex.CondTest.CondOff= sample2time(ex.(daqType).digital(1).time(2:2:end),ex.(daqType).fs,dataset.secondperunit)-ex.SufICI;
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

%% Align Scanbox frames with spike2 digital data. Peichao
    function [ex] = parsescanbox(ex)
        
        sbx = dataset.sbx.info;
        sbxFrameId = sbx.frame;
        sbxLineId = sbx.line;
        sbxLineNum = sbx.sz(1);
        sbxFrameId = sbxFrameId + sbxLineId./sbxLineNum;
        frameinTrial =  diff(sbxFrameId);
        frameinTrial = frameinTrial(1:2:end);
        timeinTrial = ex.CondTest.TrialOff-ex.CondTest.TrialOn;  % sec per trial 
        ex.timeperFrame = mean(timeinTrial'./frameinTrial);  % sec per frame
        ex.trialOnTimeCorrect =  ex.CondTest.TrialOn + (0.5-(ex.CondTest.TrialOn - floor(ex.CondTest.TrialOn))) * ex.timeperFrame;  % corrected on time
        ex.trialOffTimeCorrect =  ex.CondTest.TrialOff + (0.5-(ex.CondTest.TrialOff - floor(ex.CondTest.TrialOff))) * ex.timeperFrame; % corrected off time
        ex.firstScan = ex.CondTest.TrialOn(1) - sbxFrameId(1) * ex.timeperFrame + 0.2;  % time of start first frame.
        ex.frameTimeSer = [1:sbx.totalFrame].* ex.timeperFrame + ex.firstScan - ex.timeperFrame/2;      % time series for every frame in whole recording.
        ex.timeShiftOn = ex.frameTimeSer(sbx.frame(1:2:end)) - ex.trialOnTimeCorrect;
        ex.timeShiftOff = ex.frameTimeSer(sbx.frame(2:2:end)) - ex.trialOffTimeCorrect;
    end
end