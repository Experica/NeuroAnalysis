function ex = parseex( ex )
%PARSEEX Reads a log matrix from visstim and creates a new
%parameters structure

import NeuroAnalysis.VisStim.*

% Experimental parameters
[ex.Subject_ID, ex.File_ID, ex.RecordSite, ex.ID] = parseFilePath(ex.source);
ex.RecordSession = '';
ex.EnvParam = [];

% Monitor resolution, size, etc.
if isfield(ex.raw, 'MonRes')
    ex.EnvParam.MonitorResolution = reshape(sscanf(ex.raw.MonRes,'%f x %f'), 1, 2);
else
    ex.EnvParam.MonitorResolution = [NaN NaN];
end

if isfield(ex.raw, 'MonitorDistance')
    ex.EnvParam.ScreenToEye = ex.raw.MonitorDistance;
else
    ex.EnvParam.ScreenToEye = NaN;
end

if isfield(ex.raw, 'MonSize') && ischar(ex.raw.MonSize)
    ex.EnvParam.MonitorDegrees = reshape(sscanf(ex.raw.MonSize,'%f x %f'), 1, 2);
    ppd = mean(ex.EnvParam.MonitorResolution ./ ex.EnvParam.MonitorDegrees);
    ppcm = ppd/(2*ex.EnvParam.ScreenToEye*tand(0.5));
    ex.EnvParam.MonitorDiagonal = round(sqrt((...
        ex.EnvParam.MonitorResolution(1)^2 + ...
        ex.EnvParam.MonitorResolution(2)^2))/ppcm);
elseif isfield(ex.raw, 'MonSize')
    ex.EnvParam.MonitorDiagonal = ex.raw.MonSize*2.54; % convert to cm
else
    ex.EnvParam.MonitorDiagonal = NaN;
end

if isfield(ex.raw, 'Latency')
    ex.Latency = ex.raw.Latency;
else
    % Determine average latency of the monitor used based on the experiment
    if strcmp(ex.Subject_ID(1),'R')
        unitNo = sscanf(ex.RecordSite, 'Unit%d');
        if strcmp(ex.Subject_ID,'R1513') && unitNo <= 3
            ex.Latency = 0.024; % Samsung LCD
        elseif strcmp(ex.Subject_ID, 'R1504') || strcmp(ex.Subject_ID,'R1506') || ...
                strcmp(ex.Subject_ID,'R1508') || strcmp(ex.Subject_ID,'R1510') ...
                || strcmp(ex.Subject_ID,'R1511') || strcmp(ex.Subject_ID,'R1512')
            ex.Latency = 0.024;
        else
            ex.Latency = 0.141; % RCA tv
        end
    else
        ex.Latency = 0; % CRT
    end
end

% RF position
if isfield(ex.raw, 'RFposition')
    Position = ex.raw.RFposition;
elseif isfield(ex.raw, 'Xposition') && isfield(ex.raw, 'Yposition')
    Position = [ex.raw.Xposition,ex.raw.Yposition];
else
    Position = [0,0];
end
ex.EnvParam.Position = Position;

% BG Color
if isfield(ex.raw, 'blanks')
    ex.EnvParam.BGColor = [repmat(ex.raw.blanks/255, 1, 3) 1];
else
    ex.EnvParam.BGColor = [0.5 0.5 0.5 1];
end

% Grating type
if isfield(ex.raw, 'MakeSquare') && ~contains(ex.ID, 'Bar') && ...
        ~contains(ex.ID, 'Circle') && ~contains(ex.ID, 'Latency') && ...
        ~strcmp(ex.ID, 'RFmap') && ~strcmp(ex.ID, 'WholeScreenMapLP')
    if ex.raw.MakeSquare
        ex.EnvParam.GratingType = 'Square';
    else
        ex.EnvParam.GratingType = 'Sinusoidal';
    end
end

% Name changes
if isfield(ex.raw, 'SF') && ~contains(ex.ID, 'Spatial')
    ex.EnvParam.SpatialFreq = ex.raw.SF;
end
if isfield(ex.raw, 'TF') && ~contains(ex.ID, 'Temporal')
    ex.EnvParam.TemporalFreq = ex.raw.TF;
end
if isfield(ex.raw, 'StimSize') && ~contains(ex.ID, 'Aperture') && ...
        ~contains(ex.ID, 'Apt')
    ex.EnvParam.Diameter = ex.raw.StimSize;
    ex.EnvParam.Size = [ex.raw.StimSize, ex.raw.StimSize];
end
if isfield(ex.raw, 'XbarSize') && isfield(ex.raw, 'YbarSize') && ...
        contains(ex.ID, 'Bar')
    ex.EnvParam.Size = [ex.raw.XbarSize, ex.raw.YbarSize];
end
if isfield(ex.raw, 'C') && ~contains(ex.ID, 'Con')
    ex.EnvParam.Contrast = ex.raw.C;
end
if isfield(ex.raw, 'Ori') && ~contains(ex.ID, 'Ori') && ...
        ~contains(ex.ID, 'Pattern')
    ex.EnvParam.Ori = ex.raw.Ori;
end
if isfield(ex.raw, 'GridSize') && contains(ex.ID, 'RF')
    ex.EnvParam.GridSize = ex.raw.GridSize;
end

% Check for empty data matrices
if ~isfield(ex.raw, 'data') || isempty(ex.raw.data)
    ex.CondTest = [];
    ex.CondRepeat = 0;
    disp(['no trials for ',ex.source]);
    return
end

ex.CondSampling = 'UniformWithoutReplacement';
ex.BlockSampling = 'UniformWithoutReplacement';
ex.CondRepeat = 0;
ex.BlockRepeat = 1;
ex.BlockParam = [];

% Stim Parameters
if isfield(ex.raw, 'StimDuration')
    ex.CondDur = ex.raw.StimDuration;
else
    ex.CondDur = NaN;
end

if isfield(ex.raw, 'StimInterval')
    ex.SufICI = ex.raw.StimInterval;
else
    ex.SufICI = NaN;
end

ex.PreICI = 0;
ex.PreITI = 0;
ex.TrialDur = 0;
ex.SufITI = 0;
ex.PreIBI = 0;
ex.BlockDur = 0;
ex.SufIBI = 0;

% Fill out Data table
data = ex.raw.data;
CondTest = [];
CondTestCond = [];
switch ex.ID
    case {'velocity', 'VelocityConstantCycles', ...
            'VelocityConstantCyclesBar'}
        CondTest.StimOn = data(:,2);
        CondTest.StimOff = data(:,2) + data(:,8);
        CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,10);
        CondTestCond.Velocity = data(:,3);
    case 'Looming'
        CondTest.StimOn = data(:,2);
        CondTest.StimOff = data(:,2) + data(:,5);
        CondTest.CondRepeat = data(:,6);
        CondTest.CondIndex = data(:,7);
        CondTestCond.Velocity = data(:,3);
    case 'LatencyTest'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,4);
        CondTest.CondRepeat = data(:,3);
        CondTest.CondIndex = ones(size(data,1),1);
        CondTestCond.On = ones(size(data,1),1);
    case 'LaserGratings'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,2);
        CondTest.CondRepeat = data(:,4);
        CondTest.CondIndex = data(:,5);
        CondTestCond.On = ones(size(data,1),1);
    case 'LaserON'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,2);
        CondTest.CondRepeat = data(:,4);
        CondTest.CondIndex = ones(size(data,1),1);
        CondTestCond.On = ones(size(data,1),1);
    case 'PatternMotion'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,8);
        CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,11);
        CondTestCond.Orientation = data(:,6);
    case 'RFmap'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + 0.05;
        %CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,2);
        CondTestCond.X_Position = data(:,3);
        CondTestCond.Y_Position = data(:,4);
    case 'CatRFdetailed'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + ex.CondDur;
        %CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,2);
        CondTestCond.X_Position = data(:,4);
        CondTestCond.Y_Position = data(:,5);
        CondTestCond.Color = data(:,3);
    case {'CatRFfast10x10'}
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,6);
        CondTest.CondRepeat = data(:,7);
        CondTest.CondIndex = data(:,2);
        CondTestCond.X_Position = data(:,4);
        CondTestCond.Y_Position = data(:,5);
        CondTestCond.Color = data(:,3);
    case {'CatRFfast'} % 'CatRFfast' is broken, condition numbers are wrong
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,6);
        CondTest.CondRepeat = data(:,7);
        % CondTest.CondIndex = data(:,2);
        CondTestCond.X_Position = data(:,4);
        CondTestCond.Y_Position = data(:,5);
        CondTestCond.Color = data(:,3);
    case {'NaturalImages', 'NaturalVideos'}
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,4);
        CondTest.CondRepeat = data(:,6);
        CondTest.CondIndex = data(:,2);
        CondTestCond.Image = data(:,2);
    case {'WholeScreenMapLP'}
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,6);
        CondTest.CondRepeat = data(:,7);
        CondTest.CondIndex = data(:,2);
        CondTestCond.X_Position = data(:,4);
        CondTestCond.Y_Position = data(:,5);
        CondTestCond.Color = data(:,3);
    otherwise

        % condition, stimTime, sf, tf, apt, ori, c, dur
        if size(data,2) < 8
            error(['Unsupported log matrix for ', ex.source]);
        end
        
        CondTest.StimOn = data(:,2);
        CondTest.StimOff = data(:,2) + data(:,8);
        
        sf = data(:,3);
        tf = data(:,4);
        aperture = data(:,5);
        ori = data(:,6);
        contrast = data(:,7);
        
        % trial no
        if size(data,2) >= 9
            ex.CondTest.CondRepeat = data(:,9);
        end
        
        % velocity
        if size(data,2) >= 10
            velocity = data(:,10);
        else
            velocity = data(:,4)./data(:,3);
        end
        
        % condition name
        if contains(ex.ID, 'Ori')
            CondTestCond.Orientation = ori;
        elseif contains(ex.ID, 'Spatial')
            CondTestCond.SpatialFreq = sf;
        elseif contains(ex.ID, 'Temporal')
            CondTestCond.TemporalFreq = tf;
        elseif contains(ex.ID, 'Contrast')
            CondTestCond.Contrast = contrast;
        elseif contains(ex.ID, 'Aperture')
            CondTestCond.Diameter = aperture;
        elseif contains(ex.ID, 'Velocity')
            CondTestCond.Velocity = velocity;
        else
            error(['Unsupported log matrix for ', ex.source]);
        end
end

if isempty(CondTestCond)
    error(['No conditions for ', ex.source]);
end

% Generate new condition numbers - old ones are sometimes wrong
[Conditions, ~, conditionNo] = unique(struct2table(CondTestCond));
CondTest.CondIndex = conditionNo;
ex.Cond = table2struct(Conditions, 'ToScalar', true);

ex.CondTest = CondTest;
ex.CondTestCond = CondTestCond;

% Trim CondTest
ctnames = fieldnames(ex.CondTest);
nct = min(cellfun(@(x)length(ex.CondTest.(x)),ctnames));
for i=1:length(ctnames)
    ex.CondTest.(ctnames{i}) = ex.CondTest.(ctnames{i})(1:nct);
end

% Generate condition repeat numbers if missing
if ~isfield(ex.CondTest, 'CondRepeat')
    repeats = ones(size(unique(ex.CondTest.CondIndex),1),1);
    for i = 1:length(ex.CondTest.CondIndex)
        condIndex = ex.CondTest.CondIndex(i);
        ex.CondTest.CondRepeat(i) = repeats(condIndex);
        repeats(ex.CondTest.CondIndex(i)) = repeats(ex.CondTest.CondIndex(i)) + 1;
    end
    ex.CondRepeat = min(repeats) - 1;
else
    ex.CondRepeat = max(CondTest.CondRepeat);
end

% Convert CondIndex to integers
ex.CondTest.CondIndex = int32(ex.CondTest.CondIndex);
ex.CondTest.CondRepeat = int32(ex.CondTest.CondRepeat);

end
