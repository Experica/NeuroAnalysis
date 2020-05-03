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
[filedir,filename,ext] = fileparts(filepath);

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

metafilenames=arrayfun(@(x)x.name,dir(fullfile(filedir,[filename,'*.meta'])),'uniformoutput',0); % Find meta file from the name of analyzer file
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
    for i=1:length(metafilenames)
        m = metafilenames{i};
        mfile = fullfile(filedir,m);
        if contains(m,'imec.ap')
            dataset.ap.metafile = mfile;
            dataset.ap.meta = ReadMeta(mfile);
        elseif contains(m,'imec.lf')
            dataset.lf.metafile = mfile;
            dataset.lf.meta = ReadMeta(mfile);
        elseif contains(m,'nidq')
            dataset.nidq.metafile = mfile;
            dataset.nidq.meta = ReadMeta(mfile);
        end
    end
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
    for dss={'ap','lf','nidq'}
        ds=dss{:};
        if isfield(dataset,ds)
            meta = dataset.(ds).meta;
            meta.from = ds;
            
            samefolderbin = strrep(dataset.(ds).metafile,'meta','bin');
            if(exist(samefolderbin,'file')==2)
                meta.fileName = samefolderbin;
            elseif(exist(meta.fileName,'file')~=2)
                warning([upper(ds),' Stream Binaray File: ',meta.fileName,' not found.']);
                meta.fileName = '';
            end
            
            if strcmp(meta.typeThis,'imec')
                meta.fs = meta.imSampRate;
                % factor for converting 16-bit file data to voltage
                meta.fi2v = meta.imAiRangeMax / 512;
                meta.snsApLfSy = int64(meta.snsApLfSy);
                meta.acqApLfSy = int64(meta.acqApLfSy);
                % imec probe version
                if isfield(meta, 'imDatPrb_type')
                    meta.probeversion = str2num(meta.imDatPrb_type);
                else
                    meta.probeversion = 0; % Phase3A probe
                end
                % No. of channels multiplexed into one ADC
                switch (meta.probeversion)
                    case 2
                    otherwise
                        nchmx = 12; % Phase3A - 1.0 are all 12
                end
                nch = meta.acqApLfSy(1);
                dmxgroup = cell(nchmx,1);
                for g = 1:nchmx
                    dmxgroup{g} = g + (0:nchmx:nch-nchmx);
                end
                meta.nchmx = nchmx;
                meta.dmxgroup = dmxgroup;
                % imec readout table
                C = textscan(meta.imroTbl, '(%d %d %d %d %d', ...
                    'EndOfLine', ')', 'HeaderLines', 1 );
                meta.roch = int64(cell2mat(C(1)));
                meta.robank = int64(cell2mat(C(2)));
                meta.rorefch = int64(cell2mat(C(3)));
                meta.roapgain = cell2mat(C(4));
                meta.rolfgain = cell2mat(C(5));
                % 1-based index of reference IDs for probe option
                if meta.imProbeOpt<4
                    meta.refch = int64([36, 75, 112, 151, 188, 227, 264, 303, 340, 379]+1);
                else
                    meta.refch = int64([36, 75, 112, 151, 188, 227, 264]+1);
                end
                % bad channels
                cpn = ['SN_',num2str(meta.imProbeSN)];
                if isstruct(probeinfo) && isfield(probeinfo,cpn)
                    badch = probeinfo.(cpn).(ds).badch;
                    for i = 1: length(badch)
                        if ischar(badch{i})
                            badch{i} = str2num(badch{i});
                        end
                    end
                    meta.badch = int64(unique(cell2mat(badch))+1);
                end
                % exclude channels
                if strcmp(ds,'ap')
                    if meta.rorefch(1)==0
                        excludechans = meta.refch;
                    else
                        excludechans = meta.rorefch(1);
                    end
                else
                    excludechans = meta.refch;
                end
                if isfield(meta,'badch')
                    excludechans = union(excludechans,meta.badch);
                end
                meta.excludechans = excludechans;
                % probe channel spacing[x,y,z] in um
                meta.probespacing = [32,20,0];
                if isfield(meta,'snsShankMap')
                    header = int64(str2num(regexp(meta.snsShankMap,'([1-9,]*)','match','once')));
                    meta.nshank = header(1);
                    meta.ncol= header(2);
                    meta.nrow = header(3);
                    C = textscan(meta.snsShankMap, '(%d:%d:%d:%*s', ...
                        'EndOfLine', ')', 'HeaderLines', 1 );
                    meta.savedshanks = int64(cell2mat(C(1))+1);
                    meta.savedcols = int64(cell2mat(C(2))+1);
                    meta.savedrows = int64(cell2mat(C(3))+1);
                end
            else
                meta.fs = meta.niSampRate;
                meta.fi2v = meta.niAiRangeMax / 32768;
                meta.snsMnMaXaDw = int64(meta.snsMnMaXaDw);
            end
            
            % Return original channel IDs, because the ith channel in the file isn't necessarily
            % the ith acquired channel, so it could be used to index ith saved to original.
            if ischar(meta.snsSaveChanSubset) && strcmp(meta.snsSaveChanSubset, 'all')
                meta.savedchans = int64(1:meta.nSavedChans);
            else
                meta.savedchans = int64(meta.snsSaveChanSubset+1);
            end
            
            dataset.(ds).meta = meta;
        end
    end
    if isfield(dataset,'ap') && ~isempty(dataset.ap.meta.fileName)
        nsample=double(dataset.ap.meta.nFileSamp);
        nch = double(dataset.ap.meta.nSavedChans);
        if dataset.ap.meta.snsApLfSy(3)>0
            binmap = memmapfile(dataset.ap.meta.fileName,'Format',{'uint16',[nch,nsample],'ap'});
            digital=NeuroAnalysis.Base.parsedigitalbitinanalog(binmap.Data.ap(nch,:),nsample,16);
            if ~isempty(digital)
                for i=1:length(digital)
                    digital(i).time = NeuroAnalysis.Base.sample2time(digital(i).time,dataset.ap.meta.fs,dataset.secondperunit);
                end
                dataset.ap.digital = digital;
            end
        end
        
        if iscell(batchexportcallback) && issortconcat && ~strcmp(spikesorting,'None')
            batchexportcallback{1}=['NeuroAnalysis.SpikeGLX.',spikesorting];
        else
            switch spikesorting
                case 'KiloSort'
                    dataset = NeuroAnalysis.SpikeGLX.KiloSort(dataset);
            end
        end
    end
    disp('Preparing SpikeGLX Data:    Done.');
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
