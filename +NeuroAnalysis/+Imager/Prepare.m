function [ dataset ] = Prepare( filepath,varargin )
%PREPARE Read Imager meta file and prepare dataset
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'filepath');
addParameter(p,'ProgressBar',1,@(x)isa(x,'logical'));
addParameter(p,'IsReadAll',0,@(x)isa(x,'logical'));
addParameter(p,'exportdir','')
parse(p,filepath,varargin{:});
filepath = p.Results.filepath;
isprogressbar = p.Results.ProgressBar;
isreadall = p.Results.IsReadAll;
exportdir = p.Results.exportdir;
%% Prepare all data files
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

metafilepath =fullfile(filedir,[filename, '.meta']);
if(exist(metafilepath,'file')==2)
    isimagerdata = true;
else
    warning('No Imager Meta File:    %s', filename);
    return;
end

%% Read Imager Meta Data
if isimagerdata
    disp(['Reading Imager Meta File:    ',filename,'.meta    ...']);
    dataset = struct;
    dataset.meta = yaml.ReadYaml(metafilepath);
    dataset.source = filename;
    dataset.secondperunit = secondperunit;
    dataset.sourceformat = 'Imager';
    if ~isempty(exportdir)
        dataset.filepath = fullfile(exportdir,[filename,'.mat']);
    end
    disp(['Reading Imager Meta File:    ',filename,'.meta    Done.']);
end

%% Organize Data Files
fmt = dataset.meta.DataFormat;
width = dataset.meta.ImageFormat.Width;
height = dataset.meta.ImageFormat.Height;
pixfmt = dataset.meta.ImageFormat.PixelFormat;
file = arrayfun(@(x)x.name,dir(fullfile(filedir,[filename,'*.',fmt])),'uniformoutput',0);
ts = cellfun(@(x)regexp(x,[filename,'-Epoch(\d*)-Frame(\d*)[.]',fmt],'tokens','once'),file,'uniformoutput',0);
epoch = cellfun(@(x)str2double(x{1}),ts);
frame = cellfun(@(x)str2double(x{2}),ts);
file = cellfun(@(x)fullfile(filedir,x),file,'uniformoutput',0);

isnaturalorder=true;
checkorder = @(x) all(arrayfun(@(i)i==1,diff(unique(x,'sorted'))));
if ~checkorder(epoch)
    warning([filename, ':    Epoch Not Incremental by 1.']);
    isnaturalorder=false;
end
if ~checkorder(frame)
    warning([filename, ':    Frame Not Incremental by 1.']);
    isnaturalorder=false;
end
dataset.epochframeinorder = isnaturalorder;

ftable = table(file,epoch,frame);
ftable = sortrows(ftable,'frame');
nfile = size(ftable,1);

if isnaturalorder
    ue = unique(epoch,'sorted');
    es = cell(length(ue),1);
    for i=1:length(ue)
        es{i} = ftable.file(ftable.epoch==ue(i));
    end
    dataset.imagefile = es;
else
    dataset.imagefile = table2struct(ftable);
end

%% Read Data
if isreadall && isnaturalorder
    if isprogressbar
        hwaitbar = waitbar(0,'Reading Image Files ...');
    end
    ne = length(dataset.imagefile);
    dataset.image = cell(ne,1);
    c=0;
    for i = 1:ne
        nf = length(dataset.imagefile{i});
        dataset.image{i} = zeros(height,width,nf);
        for j = 1:nf
            dataset.image{i}(:,:,j) = NeuroAnalysis.Imager.readraw(dataset.imagefile{i}{j},width,height,pixfmt);
            c=c+1;
            if isprogressbar
                waitbar(c/nfile,hwaitbar);
            end
        end
    end
    if isprogressbar
        close(hwaitbar);
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
    % Imager data is useless without experimental data
    if ~isfield(dataset, 'ex') || isempty(dataset.ex)
        dataset = [];
    end
end

end