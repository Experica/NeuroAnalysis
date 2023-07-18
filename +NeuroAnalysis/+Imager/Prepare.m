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
dataset.meta.ImageFormat.Width = int64(dataset.meta.ImageFormat.Width);
dataset.meta.ImageFormat.Height = int64(dataset.meta.ImageFormat.Height);
pixfmt = dataset.meta.ImageFormat.PixelFormat;

files = dir(filedir);
diridx = arrayfun(@(x)x.isdir && ~strcmp(x.name,'.') && ~strcmp(x.name,'..'),files);
if any(diridx) % epochs saved in subfolders
    epochnames = arrayfun(@(x)x.name,files(diridx),'uniformoutput',0);
    ts = cellfun(@(x)regexp(x,'Epoch(\d*)','tokens','once'),epochnames,'uniformoutput',0);
    epochs = cellfun(@(x)str2double(x{1}),ts);
    [epochs,ei]=sort(epochs);
    epochnames = epochnames(ei);
    
    nepoch = length(epochnames);
    epochfiles = cell(nepoch,1);
    epochframes = cell(nepoch,1);
    for i=1:nepoch
        files = dir(fullfile(filedir,epochnames{i},[filename,'*.',fmt]));
        names = arrayfun(@(x)x.name,files,'uniformoutput',0);
        ts = cellfun(@(x)regexp(x,[filename,'-Frame(\d*)[.]',fmt],'tokens','once'),names,'uniformoutput',0);
        frames = cellfun(@(x)str2double(x{1}),ts);
        paths = cellfun(@(x)fullfile(filedir,epochnames{i},x),names,'uniformoutput',0);
        [frames,fi] = sort(frames);
        epochframes{i}=frames;
        epochfiles{i}=paths(fi);
    end
    
    isnaturalorder=true;
    if ~all(arrayfun(@(i)i==1,diff(epochs)))
        warning([filename, ':    Epoch Not Incremental by 1']);
        isnaturalorder=false;
    end
    for i=1:nepoch
        if ~all(arrayfun(@(i)i==1,diff(epochframes{i})))
            warning([filename, ':    Frame Not Incremental by 1 in Epoch ', num2str(epochs(i))]);
            isnaturalorder=false;
        end
    end
    dataset.epochframeinorder = isnaturalorder;
    
    nfile = sum(cellfun(@(x)length(x),epochframes));
    dataset.imagenepoch = int64(nepoch);
    if isnaturalorder
        dataset.imagefile = epochfiles;
    else
        dataset.imagefile.epoch = epochnames;
        dataset.imagefile.epochframe = epochframes;
        dataset.imagefile.epochfile = epochfiles;
    end
else % epochs are flatten saving
    file = arrayfun(@(x)x.name,dir(fullfile(filedir,[filename,'*.',fmt])),'uniformoutput',0);
    ts = cellfun(@(x)regexp(x,[filename,'-Epoch(\d*)-Frame(\d*)[.]',fmt],'tokens','once'),file,'uniformoutput',0);
    epochnames = cellfun(@(x)str2double(x{1}),ts);
    frame = cellfun(@(x)str2double(x{2}),ts);
    file = cellfun(@(x)fullfile(filedir,x),file,'uniformoutput',0);
    
    isnaturalorder=true;
    checkorder = @(x) all(arrayfun(@(i)i==1,diff(unique(x,'sorted'))));
    if ~checkorder(epochnames)
        warning([filename, ':    Epoch Not Incremental by 1.']);
        isnaturalorder=false;
    end
    if ~checkorder(frame)
        warning([filename, ':    Frame Not Incremental by 1.']);
        isnaturalorder=false;
    end
    dataset.epochframeinorder = isnaturalorder;
    
    ftable = table(file,epochnames,frame);
    ftable = sortrows(ftable,'frame');
    nfile = size(ftable,1);
    
    if isnaturalorder
        ue = unique(epochnames,'sorted');
        ne = length(ue);
        es = cell(ne,1);
        for i=1:ne
            es{i} = ftable.file(ftable.epoch==ue(i));
        end
        dataset.imagenepoch = int64(ne);
        dataset.imagefile = es;
    else
        dataset.imagefile = table2struct(ftable);
    end
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