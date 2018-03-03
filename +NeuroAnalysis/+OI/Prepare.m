function [ dataset ] = Prepare( filepath,varargin )
%PREPARE Read Optical Imaging block data file and prepare dataset
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'filepath');
addOptional(p,'progress',1,@(x)isa(x,'logical'));
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
progress = p.Results.progress;
%% Open block
dataset = [];
[hfile] = fopen(filepath,'r');
if hfile == -1
    warning(['Can not open file: ',filepath]);
    return;
end
%% Parse block header
disp(['Reading OI Block File:    ',filepath,'    ...']);
% Data Integrity
dataset.imagehead.filesize = fread(hfile,1,'*long');
dataset.imagehead.checksumheader = fread(hfile,1,'*long'); % Beginning with imagehead.lenheader
dataset.imagehead.checksumdata = fread(hfile,1,'*long');

% Common to all data files
dataset.imagehead.lenheader = fread(hfile,1,'*long');
dataset.imagehead.versionid = fread(hfile,1,'*long');
dataset.imagehead.filetype = fread(hfile,1,'*long'); % RAWBLOCK_FILE(11), DCBLOCK_FILE(12), SUM_FILE(13), IMAGE_FILE(14)
dataset.imagehead.filesubtype = fread(hfile,1,'*long'); % FROM_VDAQ(11), FROM_ORA(12)
dataset.imagehead.datatype = fread(hfile,1,'*long'); % DAT_UCHAR(11), DAT_USHORT(12), DAT_LONG(13), DAT_FLOAT(14)
dataset.imagehead.sizeof = fread(hfile,1,'*long'); % e.g. sizeof(long), sizeof(float)
dataset.imagehead.framewidth = fread(hfile,1,'*long');
dataset.imagehead.frameheight = fread(hfile,1,'*long');
dataset.imagehead.nframesperstim = fread(hfile,1,'*long'); % data frames
dataset.imagehead.nstimuli = fread(hfile,1,'*long');
dataset.imagehead.initialxbinfactor = fread(hfile,1,'*long'); % from data acquisition
dataset.imagehead.initialybinfactor = fread(hfile,1,'*long'); % from data acquisition
dataset.imagehead.xbinfactor = fread(hfile,1,'*long'); % this file
dataset.imagehead.ybinfactor = fread(hfile,1,'*long'); % this file
dataset.imagehead.username = fread(hfile,32,'*char');
dataset.imagehead.recordingdate = fread(hfile,16,'*char');
dataset.imagehead.x1roi = fread(hfile,1,'*long');
dataset.imagehead.y1roi = fread(hfile,1,'*long');
dataset.imagehead.x2roi = fread(hfile,1,'*long');
dataset.imagehead.y2roi = fread(hfile,1,'*long');

% Locate data and ref frames
dataset.imagehead.stimoffs = fread(hfile,1,'*long');
dataset.imagehead.stimsize = fread(hfile,1,'*long');
dataset.imagehead.frameoffs = fread(hfile,1,'*long');
dataset.imagehead.framesize = fread(hfile,1,'*long');
dataset.imagehead.refoffs = fread(hfile,1,'*long'); % Imager 3001 has no ref
dataset.imagehead.refsize = fread(hfile,1,'*long'); % these fields will be 0
dataset.imagehead.refwidth = fread(hfile,1,'*long');
dataset.imagehead.refheight = fread(hfile,1,'*long');

% Common to data files that have undergone some form of
% "compression" or "summing"; i.e. The data in the current
% file may be the result of having summed blocks 'a'-'f', frames 1-7
dataset.imagehead.whichblocks = fread(hfile,16,'*ushort'); % 256 bits => max of 256 blocks per experiment
dataset.imagehead.whichframes = fread(hfile,16,'*ushort'); % 256 bits => max of 256 frames per condition

% Data analysis
dataset.imagehead.loclip = fread(hfile,1,'float');
dataset.imagehead.hiclip = fread(hfile,1,'float');
dataset.imagehead.lopass = fread(hfile,1,'*long');
dataset.imagehead.hipass = fread(hfile,1,'*long');
dataset.imagehead.operationsperformed = fread(hfile,64,'*char');

% Ora specific - not needed by Vdaq
dataset.imagehead.magnification = fread(hfile,1,'float');
dataset.imagehead.gain = fread(hfile,1,'*ushort');
dataset.imagehead.wavelength = fread(hfile,1,'*ushort');
dataset.imagehead.exposuretime = fread(hfile,1,'*long');
dataset.imagehead.nrepetitions = fread(hfile,1,'*long'); % number of repetitions
dataset.imagehead.acquisitiondelay = fread(hfile,1,'*long'); % delay of DAQ relative to Stim-Go
dataset.imagehead.interstiminterval = fread(hfile,1,'*long'); % time interval between Stim-Go's
dataset.imagehead.creationdate = fread(hfile,16,'*char');
dataset.imagehead.datafilename = fread(hfile,64,'*char');
dataset.imagehead.orareserved = fread(hfile,256,'*char');

% Vdaq-specific
dataset.imagehead.includesrefframe = fread(hfile,1,'*long'); % 0 or 1
dataset.imagehead.listofstimuli = fread(hfile,256,'*char');
dataset.imagehead.nframesperdataframe = fread(hfile,1,'*long');
dataset.imagehead.ntrials = fread(hfile,1,'*long');
dataset.imagehead.scalefactor = fread(hfile,1,'*long'); % NFramesAvgd * Bin * Trials
dataset.imagehead.meanampgain = fread(hfile,1,'float');
dataset.imagehead.meanampdc = fread(hfile,1,'float');
dataset.imagehead.begbaselineframeno = fread(hfile,1,'*int8'); % SUM-FR/DC File (i.e. compressed)
dataset.imagehead.endbaselineframeno = fread(hfile,1,'*int8'); % SUM-FR/DC File (i.e. compressed)
dataset.imagehead.begactivityframeno = fread(hfile,1,'*int8'); % SUM-FR/DC File (i.e. compressed)
dataset.imagehead.endactivityframeno = fread(hfile,1,'*int8'); % SUM-FR/DC File (i.e. compressed)
dataset.imagehead.digitizerbits = fread(hfile,1,'*int8'); % cam_GetGrabberBits
dataset.imagehead.activesystemid = fread(hfile,1,'*int8'); % core_ActiveSystemID()
dataset.imagehead.dummy2 = fread(hfile,1,'*int8');
dataset.imagehead.dummy3 = fread(hfile,1,'*int8');
dataset.imagehead.x1superpix = fread(hfile,1,'*long');
dataset.imagehead.y1superpix = fread(hfile,1,'*long');
dataset.imagehead.x2superpix = fread(hfile,1,'*long');
dataset.imagehead.y2superpix = fread(hfile,1,'*long');
dataset.imagehead.frameduration = fread(hfile,1,'float');
dataset.imagehead.validframes = fread(hfile,1,'*long');
dataset.imagehead.vdaqreserved = fread(hfile,224,'*char');

% SYSTEMTIME is 8 WORDS (16 bytes)
dataset.imagehead.timeblockstart = fread(hfile,8,'*int16');
dataset.imagehead.timeblockend = fread(hfile,8,'*int16');
dataset.imagehead.user = fread(hfile,224,'*char');
dataset.imagehead.comment = fread(hfile,256,'*char');
%% Prepare block header
dataset.imagehead.listofstimuli = int32(str2double(strsplit(dataset.imagehead.listofstimuli')));
%% Prepare experimental data
[~,fname,~]=fileparts(filepath);
ts = regexp(fname,'^([A-Za-z0-9]+)_(E[0-9]+)B[0-9]+','tokens','once');
dataset.Subject_ID=ts{1};
dataset.RecordSession=ts{2};
%% Read block
switch dataset.imagehead.datatype
    case 11 % DAT_UCHAR
        datatype = '*uchar';
    case 12 % DAT_USHORT
        datatype = '*ushort';
    case 13 % DAT_LONG
        datatype = '*long';
    case 14 % DAT_FLOAT
        datatype = '*single';
end
if progress
    hwaitbar = waitbar(0,'Reading OI Block File ...');
end
dataset.image = cell(dataset.imagehead.ntrials,dataset.imagehead.nstimuli);
for i = 1:dataset.imagehead.ntrials
    for s = 1:dataset.imagehead.nstimuli
        dataset.image{i,s} = zeros(dataset.imagehead.frameheight,dataset.imagehead.framewidth,dataset.imagehead.nframesperstim);
        for d = 1:dataset.imagehead.nframesperstim
            ti = fread(hfile,[dataset.imagehead.framewidth,dataset.imagehead.frameheight],datatype);
            dataset.image{i,s}(:,:,d) = ti';
            if progress
                waitbar((i*s*d)/(dataset.imagehead.ntrials*dataset.imagehead.nstimuli*dataset.imagehead.nframesperstim),hwaitbar);
            end
        end
    end
end
if progress
    close(hwaitbar);
end
if ~isempty(dataset)
    dataset.source = filepath;
    dataset.sourceformat = 'OI';
end
if (ftell(hfile)~=dataset.imagehead.filesize)
    warning('Block Is Not Completely Read Out.');
end
fclose(hfile);
disp('Reading File:    Done.');
end