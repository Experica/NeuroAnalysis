function [ dataset ] = Prepare( filepath,varargin )
%PREPARE Read SpikeGLX meta data, experimental data and prepare dataset
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'filepath');
addParameter(p,'SpikeSorting','None')
addParameter(p,'IsConcat',0)
addParameter(p,'exportdir','')
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
spikesorting = p.Results.SpikeSorting;
issortconcat = p.Results.IsConcat;
exportdir = p.Results.exportdir;
global batchexportcallback
%% Prepare all data files
probeinfo=[];
[thisdir,~,~] = fileparts(mfilename('fullpath'));
probefilepath = fullfile(thisdir,'NeuropixelsProbeInfo.yaml');
if(exist(probefilepath,'file')==2)
    probeinfo = yaml.ReadYaml(probefilepath);
end

dataset = [];
secondperunit=1;
[filedir,filename,~] = fileparts(filepath);

isexdata = false;
exfilepath =fullfile(filedir,[filename '.yaml']);
if(exist(exfilepath,'file')==2)
    isexdata = true;
    secondperunit=0.001;
end

isstimulatordata = false;
analyzerfilepath =fullfile(filedir,[filename, '.analyzer']);
if(exist(analyzerfilepath,'file')==2)
    isstimulatordata = true;
end

metafilenames=arrayfun(@(x)x.name,dir(fullfile(filedir,[filename,'*.meta'])),'uniformoutput',0);
if ~isempty(metafilenames)
    isspikeglxdata = true;
else
    warning('No SpikeGLX Meta Files:    %s', filename);
    return;
end

%% Read SpikeGLX Meta Data
if(isspikeglxdata)
    disp(['Reading SpikeGLX Meta Files:    ',filename,'*.meta    ...']);
    dataset=struct;
    imecindex={};
    for i=1:length(metafilenames)
        m = metafilenames{i};
        mfile = fullfile(filedir,m);
        ts = regexp(m,'\w*[.](imec|nidq)(\d*)[.]?(ap|lf)?[.]meta','tokens','once');
        f = ts{1};
        if strcmp(f,'imec')
            imecindex=[imecindex,ts(2)];
            f = [ts{3},ts{2}];
        end
        dataset.(f).metafile = mfile;
        dataset.(f).meta = NeuroAnalysis.SpikeGLX.readmeta(mfile);
    end
    dataset.imecindex = unique(imecindex);
    disp(['Reading SpikeGLX Meta Files:    ',filename,'*.meta    Done.']);
end

if ~isempty(dataset)
    dataset.source = filename;
    dataset.secondperunit = secondperunit;
    dataset.sourceformat = 'SpikeGLX';
    if ~isempty(exportdir)
        dataset.filepath = fullfile(exportdir,[filename,'.mat']);
    end
end
%% Prepare SpikeGLX data
if ~isempty(dataset)
    disp('Preparing SpikeGLX Data:    ...');
    fs = fields(dataset);
    fs = fs(cellfun(@(x)~isempty(regexp(x,'(nidq|ap|lf)(\d*)', 'once')),fs));
    for j=1:length(fs)
        d=fs{j};
        meta = dataset.(d).meta;
        meta.from = d;
        
        samefolderbin = strrep(dataset.(d).metafile,'meta','bin');
        if(exist(samefolderbin,'file')==2)
            meta.fileName = samefolderbin;
        elseif(exist(meta.fileName,'file')~=2)
            warning([upper(d),' Stream Binaray File: ',meta.fileName,' not found.']);
            meta.fileName = '';
        end
        dataset.(d).meta = NeuroAnalysis.SpikeGLX.parsemeta(meta,d,probeinfo);
        
        % parse digital data
        if ~startsWith(d,'lf') && ~isempty(meta.fileName)
            nsample=double(meta.nFileSamp);
            nch = double(meta.nSavedChans);
            binmap = memmapfile(meta.fileName,'Format',{'uint16',[nch,nsample],'d'});
            digital = NeuroAnalysis.Base.parsedigitalbitinanalog(binmap.Data.d(nch,:),nsample,16);
            if ~isempty(digital)
                for i=1:length(digital)
                    digital(i).time = NeuroAnalysis.Base.sample2time(digital(i).time,meta.fs,dataset.secondperunit);
                end
                dataset.(d).digital = digital;
            end
            
            if startsWith(d,'ap')
                if iscell(batchexportcallback) && issortconcat && ~strcmp(spikesorting,'None')
                    batchexportcallback{1}=['NeuroAnalysis.SpikeGLX.',spikesorting];
                else
                    switch spikesorting
                        case 'KiloSort'
                            dataset = NeuroAnalysis.SpikeGLX.KiloSort(dataset,d);
                    end
                end
            end
        end
    end
    
    % digital markers, and sync timing difference between imec and nidq
    if ~isfield(dataset,'digital')
        digital=[];
        if isempty(dataset.imecindex{1}) % Phase3A, no sync
            if isfield(dataset,'ap') && isfield(dataset.ap,'digital')
                digital = dataset.ap.digital;
            end
        else % PXI system, markers in nidq
            if isfield(dataset,'nidq') && isfield(dataset.nidq,'digital')
                di = arrayfun(@(x)x.channel~=dataset.nidq.meta.syncch,dataset.nidq.digital);
                digital = dataset.nidq.digital(di);
                sync = dataset.nidq.digital(~di);
                if ~isempty(sync)
                    % clean any narrow pluses in nidq sync, here only handle dirac pluses
                    % on digital low state(may be the ground level not equal or noisy).
                    [synctime,syncdata] = NeuroAnalysis.Base.cleandigitalpluse(sync.time,sync.data,0.1,true,'for nidq sync');
                    dataset.sync = synctime;
                    dataset.syncdt = mean(diff(dataset.sync));
                    
                    for i=1:length(dataset.imecindex)
                        d = ['ap',dataset.imecindex{i}];
                        s = ['syncdiff',dataset.imecindex{i}];
                        if isfield(dataset,d) && isfield(dataset.(d),'digital')
                            ti = find(arrayfun(@(x)x.channel==dataset.(d).meta.syncch,dataset.(d).digital));
                            if isempty(ti)
                                warning('Sync Data Not Found In %s, Skip SyncDiff ...',d);
                            else
                                dataset.(s) = syncdiff(dataset.sync,dataset.(d).digital(ti).time);
                            end
                        end
                    end
                end
                if ~isempty(digital)
                    % clean any narrow pluses in nidq markers, here only handle dirac pluses
                    % on digital low state(may be the ground level not equal or noisy).
                    for i=1:length(digital)
                        [digital(i).time,digital(i).data] = NeuroAnalysis.Base.cleandigitalpluse(digital(i).time,digital(i).data,0.1,true,['for nidq digital channel ',num2str(digital(i).channel)]);
                    end
                end
            end
        end
        dataset.digital = digital;
    end
    
    disp('Preparing SpikeGLX Data:    Done.');
end

    function [d] = syncdiff(rt,tt)
        d = [];
        rn = length(rt);
        tn = length(tt);
        if rn == tn
            d = tt - rt;
        else
            if abs(rn-tn)==1
                warning('Matching %i Ref Sync Pluse to %i Target Sync Pluse ...',rn,tn);
                tt1 = tt(1);
                rt1 = rt(1);
                if abs(tt1-rt1) > dataset.syncdt/2
                    if tt1>rt1
                        d = [0,tt-rt(2:end)];
                    else
                        d = tt(2:end)-rt;
                    end
                else
                    if rn>tn
                        d = [tt-rt(1:end-1),0];
                    else
                        d = tt(1:end-1)-rt;
                    end
                end
            else
                warning('Number of Sync Pluse(%i/%i) Not Matching, Skip SyncDiff ...',rn,tn);
            end
        end
    end
%% Prepare corresponding experimental data
if ~isempty(dataset)
    if(isexdata)
        exdataset = NeuroAnalysis.Experica.Prepare(exfilepath,dataset);
        if ~isempty(exdataset)
            dataset.ex = exdataset.ex;
        end
    elseif(isstimulatordata)
        stimulatordataset = NeuroAnalysis.Stimulator.Prepare(analyzerfilepath,dataset);
        if ~isempty(stimulatordataset)
            dataset.ex = stimulatordataset;
        end
    end
    % SpikeGLX data is useless without experimental data
    if ~isfield(dataset, 'ex') || isempty(dataset.ex)
        dataset = [];
    end
end

end
