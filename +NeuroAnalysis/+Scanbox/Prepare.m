function [ dataset ] = Prepare( filepath,varargin )
%PREPARE combine some data together for analysis
% Peichao: Combine Scanbox .mat info, analyzer, spkie2 data and hartely log files if have one
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'filepath');
addParameter(p,'SpikeSorting','None')
addParameter(p,'IsConcat',0)
addParameter(p,'exportdir','')
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
% spikesorting = p.Results.SpikeSorting;
% issortconcat = p.Results.IsConcat;
% exportdir = p.Results.exportdir;

%% Prepare all data files
dataset = struct;
% dataset = [];

% filepath = 'O:\AF4\2P_data\U004\metaFiles\AF4_u004_005.analyzer';  % for test

[filedir,filename,ext] = fileparts(filepath);

isexdata = false;
exfilepath =fullfile(filedir,[filename '.yaml']);
if(exist(exfilepath,'file')==2)
    isexdata = true;
end

isstimulatordata = false;
analyzerfilepath = fullfile(filedir,[filename, '.analyzer']);
if(exist(analyzerfilepath,'file')==2)
    isstimulatordata = true;
    dataset.secondperunit = 1; % in sec
end
spike2File = dir(fullfile(filedir,[filename, '*.smr']));
spike2MatFile = dir(fullfile(filedir,[filename, '*hartley*.mat']));
isHartleydata = false;
dataset.eventSource = '';
if (~isempty(spike2File)) || (~isempty(spike2MatFile))
    isHartleydata = true;
    dataset.eventSource = 'spike2';
    spike2TimeUnit = 1e-6;  %  1 micro sec per unit, feature of CED1401 
end

% Tansform spike2 .smr to matlab .mat
if isHartleydata
    if ~exist(fullfile(spike2MatFile.folder, spike2MatFile.name), 'file')
        disp(['Reading Spike2 .smr Files:    ',spike2File.name,'    ...']);
        [time, interval, ttlData, chanName]=load_smr(fullfile(spike2File.folder, spike2File.name));
        spike2Data = [];
        for ii = 1: length(chanName)
            field = chanName{ii};
            value = ttlData(:,ii);
            spike2Data.(field) = value;
        end
        spike2Data.time = time;
        spike2Data.fs = 1/(interval * spike2TimeUnit);  % us to s, and to sampling rate
        spike2Data.chanName = chanName;
        disp(['Saving Spike2 .smr Files as .mat file     ...']);
        save(fullfile(spike2File.folder, [spike2File.name(1:end-4), '.mat']), 'spike2Data', '-v7.3');
        dataset.spike2Data = spike2Data;
        dataset.exptId = 'Hartley';
    else
        disp(['Loading Spike2 .mat Files:    ',spike2MatFile.name,'    ...']);
        aa=load(fullfile(spike2MatFile.folder, spike2MatFile.name), '-mat');
        % The following commented code is for spkie2 data saved as .mat file, not .smr file
%         chanName=fieldnames(aa); 
%         spike2Data = [];
%         for ii = 1: length(chanName)
%             field = aa.(chanName{ii}).title;
%             value = aa.(chanName{ii}).values;
%             spike2Data.(field) = value;
%         end
%         interval = aa.(chanName{1}).interval;
%         sampNum = aa.(chanName{1}).length;
%         spike2Data.time = interval*(0:(sampNum-1))';
%         spike2Data.fs = 1/(interval);  % sampling rate
%         spike2Data.chanName = chanName;
%         disp(['Saving Spike2 .mat Files as .mat file     ...']);
%         save(fullfile(spike2MatFile.folder, [spike2MatFile.name(1:end-4), 'new.mat']), 'spike2Data', '-v7.3');
        dataset.spike2Data = aa.spike2Data;
        dataset.exptId = 'Hartley';
    end

else
    dataset.exptId = '';
    warning('No Spike2 .smr  Files!!!!');
end


ani = filename(1:3);
unit = filename(6:8);
expt = filename(10:12);

filenameNew = [ani, '_', unit, '_', expt];
dirinfo = dir(fullfile(filedir, [filenameNew,'*_parse.scan']));

% dirinfo = dir(fullfile(filedir,'\*\*\', [filenameNew,'*_split.mat']));

% for i=1:length(dirinfo)
%     if length(dirinfo(i,:).name)>22
%         dirinfo(i,:)=[];
%     end
% end 

metafilenames=arrayfun(@(x)x.name,dirinfo,'uniformoutput',0); % Find meta file from the name of analyzer file

if ~isempty(metafilenames)
    issbxdata = true;
else
    warning('No Scanbox Mat Files:    %s', filename);
    return;
end

%% Read Scanbox Meta Data
if(issbxdata)
    disp(['Reading Scanbox mat Files:    ',filenameNew,'_parse.scan    ...']);
    
    for i=1:length(metafilenames)
        mfile = fullfile(dirinfo(i).folder, metafilenames{i});
        dataset.sbx.metafile = mfile;
        dataset.sbx = load(mfile,'-mat');
        disp(['Reading Scanbox mat Files:    ',mfile,'    Done.']);
        if ~isempty(dataset)
            dataset.source = metafilenames{i};
            dataset.sourceformat = 'Scanbox';
            dataset.DAQformat = 'Scanbox';
%             if ~isempty(exportdir)
            dataset.filepath = fullfile(dirinfo(i).folder,[metafilenames{i}(1:end-10),'meta.mat']);
%             exportdir = dataset.filepath;
%             end
        % Prepare Scanbox data
            disp('Preparing Scanbox Data:    ...');

            meta = dataset.sbx.info;
            meta.fileName = metafilenames{i};
            if isempty(meta.otparam)
                meta.planeNum = 1;
            elseif ~isempty(meta.otparam)
                meta.planeNum = meta.otparam(3);
            end
%             if (isfield(meta, 'frame_split'))
%                 meta.frame = metafilenames{i}(end-2:end);
%             else
%                 
%             end
            disp('Preparing Scanbox Data:    Done.');
        end
        
        % Prepare corresponding experimental data
        if ~isempty(dataset)
            if(isexdata)
                exdataset = NeuroAnalysis.Experica.Prepare(exfilepath,dataset);
                if ~isempty(exdataset)
                    dataset.ex = exdataset.ex;
                end
            elseif(isstimulatordata)
                stimulatordataset = NeuroAnalysis.Stimulator.Prepare(analyzerfilepath, dataset);
                if ~isempty(stimulatordataset)
                    dataset.ex = stimulatordataset;
                end
            end
            % Scanbox data is useless without experimental data
            if ~isfield(dataset, 'ex') || isempty(dataset.ex)
                dataset = [];
            end
        end
        
        
   end
end
end
