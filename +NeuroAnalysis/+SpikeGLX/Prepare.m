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
            meta = setfield(meta, tag, v);
        end
        meta.nFileSamp = meta.fileSizeBytes/meta.nSavedChans/2;
        meta.fileDate = datenum(meta.fileCreateTime,'yyyy-mm-ddTHH:MM:SS');
    end

p = inputParser;
addRequired(p,'filepath');
addParameter(p,'SpikeSorting','None')
addParameter(p,'IsConcat',0)
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
spikesorting = p.Results.SpikeSorting;
issortconcat = p.Results.IsConcat;
global batchexportcallback
%% Prepare all data files
dataset = [];
secondperunit=1;
[filedir,filename,ext] = fileparts(filepath);

metafilenames=arrayfun(@(x)x.name,dir(fullfile(filedir,[filename,'*.meta'])),'uniformoutput',0);

if ~isempty(metafilenames)
    isspikeglxdata = true;
else
    warning('No SpikeGLX Meta Files:    %s', filename);
    return;
end

isexdata = false;
exfilepath =fullfile(filedir,[filename '.yaml']);
if(exist(exfilepath,'file')==2)
    isexdata = true;
    secondperunit=0.001;
end

isstimulatordata = false;
analyzerfilepath =fullfile(filedir,[filename '.analyzer']);
if(exist(analyzerfilepath,'file')==2)
    isstimulatordata = true;
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
end
%% Prepare SpikeGLX data
if ~isempty(dataset)
    disp('Preparing SpikeGLX Data:    ...');
    for dss={'ap','lf','nidq'}
        ds=dss{:};
        if isfield(dataset,ds)
            samefolderbin = strrep(dataset.(ds).metafile,'meta','bin');
            if(exist(samefolderbin,'file')==2)
                dataset.(ds).meta.fileName = samefolderbin;
            elseif(exist(dataset.(ds).meta.fileName,'file')~=2)
                warning([upper(ds),' Stream Binaray File: ',dataset.(ds).meta.fileName,' not found.']);
                dataset.(ds).meta.fileName = '';
            end
        end
    end
    if isfield(dataset,'ap') && ~isempty(dataset.ap.meta.fileName)
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
