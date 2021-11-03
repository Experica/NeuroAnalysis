function [ dataset ] = Prepare( filepath,varargin )
%PREPARE Read SpikeGLX meta data, experimental data and prepare dataset
%   Detailed explanation goes here

    function [meta] = ReadMeta(filepath)
        % Parse ini file into cell entries C{1}{i} = C{2}{i}
        fid = fopen(filepath, 'r');
        C = textscan(fid, '%[^=] = %[^\r\n]');
        fclose(fid);
        
        meta = struct();
        % Convert each cell entry into a struct entry
        for tagi = 1:length(C{1})
            tag = C{1}{tagi};
            if tag(1) == '~'
                % remake tag excluding first character
                tag = sprintf('%s', tag(2:end));
            end
            v = C{2}{tagi};
            [t,r]=str2num(v);
            if r
                v=t;
            end
            meta.(tag) = v;
        end
        t = meta.fileSizeBytes/meta.nSavedChans;
        if mod(t,2) ~= 0
            warning('Binary file does not have same number of samples for all channels. Use the maximum same length for all channels.');
            meta.nFileSamp = int64(floor(t/2));
        else
            meta.nFileSamp = int64(t/2);
        end
        meta.nSavedChans = int64(meta.nSavedChans);
        meta.fileDate = datenum(meta.fileCreateTime,'yyyy-mm-ddTHH:MM:SS');
    end

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
        dataset.(f).meta = ReadMeta(mfile);
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
        
        if strcmp(meta.typeThis,'imec')
            meta.fs = meta.imSampRate;
            meta.snsApLfSy = int64(meta.snsApLfSy);
            meta.acqApLfSy = int64(meta.acqApLfSy);
            % imec probe version
            if isfield(meta, 'imDatPrb_type')
                meta.probeversion = meta.imDatPrb_type;
            else
                meta.probeversion = -1; % Phase3A
            end
            switch (meta.probeversion)
                case {21,24}
                    % factor for converting 16-bit data to voltage
                    fi2v = meta.imAiRangeMax / 8192;
                otherwise % Phase3A - 1.0
                    % probe channel spacing[x,y,z] in um
                    spacing = [32,20,0];
                    % factor for converting 16-bit data to voltage
                    fi2v = meta.imAiRangeMax / 512;
                    % No. of channels multiplexed into one ADC
                    nchmx = 12;
                    % ADC1 {0,2,4,6,8,10,12,14,16,18,20,22}
                    % ADC2 {1,3,5,7,9,11,13,15,17,19,21,23}
                    % ADC3 {24,26,28,30,32,34,36,38,40,42,44,46}
                    % ADC4 {25,27,29,31,33,35,37,39,41,43,45,47}
                    % ...
                    nch = meta.acqApLfSy(1);
                    dmxgroup = cell(nchmx,1);
                    hc = 1:2*nchmx:nch;
                    dmxgroup{1} = [hc,hc+1];
                    for g = 2:nchmx
                        dmxgroup{g} = dmxgroup{g-1}+2;
                    end
                    % reference IDs
                    if isfield(meta, 'imProbeOpt') % Phase3A
                        switch  meta.imProbeOpt
                            case 4
                                refch = int64([36, 75, 112, 151, 188, 227, 264]+1);
                            otherwise
                                refch = int64([36, 75, 112, 151, 188, 227, 264, 303, 340, 379]+1);
                        end
                        probesn = meta.imProbeSN;
                        rofmt = '(%d %d %d %d %d';
                    else % Phase3B - 1.0
                        refch = int64(191+1); % 192, 576, 960 for bank 0, 1, 2
                        probesn = meta.imDatPrb_sn;
                        rofmt = '(%d %d %d %d %d %d';
                    end
                    % imec readout table
                    C = textscan(meta.imroTbl, rofmt, ...
                        'EndOfLine', ')', 'HeaderLines', 1 );
                    meta.roch = int64(cell2mat(C(1)));
                    meta.robank = int64(cell2mat(C(2)));
                    meta.rorefch = int64(cell2mat(C(3)));
                    meta.roapgain = cell2mat(C(4));
                    meta.rolfgain = cell2mat(C(5));
            end
            meta.nchmx = nchmx;
            meta.dmxgroup = dmxgroup;
            meta.probespacing = spacing;
            meta.refch = refch;
            meta.fi2v = fi2v;
            meta.probesn = probesn;
            % bad channels
            psn = ['SN_',num2str(meta.probesn)];
            if isstruct(probeinfo) && isfield(probeinfo,psn)
                badch = probeinfo.(psn).(d).badch;
                for i = 1:length(badch)
                    if ischar(badch{i})
                        badch{i} = str2num(badch{i});
                    end
                end
                meta.badch = int64(unique(cell2mat(badch))+1);
            end
            % exclude channels, here only handle external referencing
            if all(meta.rorefch==0)
                excludechans = meta.refch;
            end
            if isfield(meta,'badch')
                excludechans = union(excludechans,meta.badch);
            end
            meta.excludechans = excludechans;
            % imec shank map for saved channels
            if isfield(meta,'snsShankMap')
                header = int64(str2num(regexp(meta.snsShankMap,'([0-9,]*)','match','once')));
                meta.nshank = header(1);
                meta.ncol= header(2);
                meta.nrow = header(3);
                C = textscan(meta.snsShankMap, '(%d:%d:%d:%*s', ...
                    'EndOfLine', ')', 'HeaderLines', 1 );
                meta.savedshanks = int64(cell2mat(C(1))+1);
                meta.savedcols = int64(cell2mat(C(2))+1);
                meta.savedrows = int64(cell2mat(C(3))+1);
                meta.nshanksaved = int64(length(unique(meta.savedshanks)));
                meta.ncolsaved = int64(length(unique(meta.savedcols)));
                meta.nrowsaved = int64(length(unique(meta.savedrows)));
            end
            meta.syncch = 6+1; % fixed bit 6 for imec sync in PXI system
        else % nidq
            meta.fs = meta.niSampRate;
            % factor for converting 16-bit data to voltage
            meta.fi2v = meta.niAiRangeMax / 32768;
            meta.snsMnMaXaDw = int64(meta.snsMnMaXaDw);
            meta.syncch = meta.syncNiChan + 1;
        end
        % Return original channel IDs, because the ith channel in the file isn't necessarily
        % the ith acquired channel, so it could be used to index ith saved to original.
        if ischar(meta.snsSaveChanSubset) && strcmp(meta.snsSaveChanSubset, 'all')
            meta.savedchans = int64(1:meta.nSavedChans);
        else
            meta.savedchans = int64(meta.snsSaveChanSubset+1);
        end
        dataset.(d).meta = meta;
        
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
                    % clean any random noise in nidq sync, here only handle random pluses on digital low state.
                    synctime = sync.time;
                    ri = find(sync.data~=0);
                    fi = ri+1;
                    if fi(end)>length(sync.data)
                        fi = fi(1:end-1);
                        ri = ri(1:end-1);
                    end
                    dt = dataset.nidq.meta.syncSourcePeriod/secondperunit/2;
                    ni = find(synctime(fi)-synctime(ri) < dt/2);
                    if ~isempty(ni)
                        warning('Clean Noisy Sync Data ...');
                        synctime([ri(ni),fi(ni)])=[];
                    end
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
