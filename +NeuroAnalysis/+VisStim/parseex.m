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
    ex.EnvParam.ScreenResolution = reshape(sscanf(ex.raw.MonRes,'%f x %f'), 1, 2);
else
    ex.EnvParam.ScreenResolution = [NaN NaN];
end

if isfield(ex.raw, 'MonitorDistance')
    ex.EnvParam.ScreenToEye = ex.raw.MonitorDistance;
else
    ex.EnvParam.ScreenToEye = NaN;
end

if isfield(ex.raw, 'MonSize') && isnumeric(ex.raw.MonSize)
    ex.EnvParam.ScreenDiagonal = ex.raw.MonSize*2.54; % convert to cm
elseif isfield(ex.raw, 'MonitorSize') && isnumeric(ex.raw.MonitorSize)
    ex.EnvParam.ScreenDiagonal = ex.raw.MonitorSize*2.54; % convert to cm
else
    ex.EnvParam.ScreenDiagonal = NaN;
end

ex.EnvParam.ScreenAspect = ex.EnvParam.ScreenResolution(1)/...
    ex.EnvParam.ScreenResolution(2);
ex.EnvParam.ScreenHalfHeight = sqrt(ex.EnvParam.ScreenDiagonal^2/...
    (1 + ex.EnvParam.ScreenAspect^2))/2;
deg = atand(ex.EnvParam.ScreenHalfHeight/ex.EnvParam.ScreenToEye);
ex.EnvParam.ScreenDegrees = [deg*ex.EnvParam.ScreenAspect*2, ...
    deg*2];

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
    ex.EnvParam.BGColor = [repmat(ex.raw.blanks, 1, 3) 1];
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
    warning(['no trials for ',ex.source]);
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
    ex.CondDur = 0;
end
if isfield(ex.raw, 'StimInterval')
    if ex.raw.StimInterval > 0.5
        ex.PreICI = 0.5;
    else
        ex.PreICI = ex.raw.StimInterval/2;
    end
    ex.SufICI = ex.raw.StimInterval - ex.PreICI;
else
    ex.PreICI = 0;
    ex.SufICI = 0;
end

ex.PreITI = 0;
ex.TrialDur = 0;
ex.SufITI = 0;
ex.PreIBI = 0;
ex.BlockDur = 0;
ex.SufIBI = 0;

% Fill out Data table
data = ex.raw.data;
CondTest = struct;
% Possible factors
Factors = struct;
switch ex.ID
    case {'velocity', 'VelocityConstantCycles', ...
            'VelocityConstantCyclesBar'}
        CondTest.StimOn = data(:,2);
        CondTest.StimOff = data(:,2) + data(:,8);
        CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,10);
        Factors.Velocity = data(:,3);
        Factors.nDots = data(:,4);
        Factors.Diameter = data(:,5);
        Factors.Ori = data(:,6);
        Factors.Contrast = data(:,7);
    case 'Looming'
        CondTest.StimOn = data(:,2);
        CondTest.StimOff = data(:,2) + data(:,5);
        CondTest.CondRepeat = data(:,6);
        CondTest.CondIndex = data(:,7);
        Factors.Velocity = data(:,3);
        Factors.Contrast = data(:,4);
    case 'LatencyTest'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,4);
        CondTest.CondRepeat = data(:,3);
        CondTest.CondIndex = ones(size(data,1),1);
        Factors.Visible = ones(size(data,1),1);
    case 'LaserGratings'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,2);
        CondTest.CondRepeat = data(:,4);
        CondTest.CondIndex = data(:,5);
        Factors.Visible = ones(size(data,1),1);
    case 'LaserON'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,2);
        CondTest.CondRepeat = data(:,4);
        CondTest.CondIndex = ones(size(data,1),1);
        Factors.Visible = ones(size(data,1),1);
    case 'PatternMotion'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,8);
        CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,11);
        Factors.Orientation = data(:,6);
        Factors.Diameter = data(:,5);
        Factors.Ori = data(:,6);
        Factors.Contrast = data(:,7);
    case 'RFmap'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + 0.05;
        %CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,2);
        Factors.Position = [data(:,3), data(:,4), zeros(size(data,1),1)];
        Factors.Diameter = data(:,5);
        ex.CondDur = 0.05;
        ex.SufICI = 0.3;
    case 'Retinotopy'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + ex.raw.StimInterval;
        %CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,2);
        Factors.Position = [data(:,3), data(:,4), zeros(size(data,1),1)];
        Factors.Size = data(:,5);
    case 'CatRFdetailed'
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + ex.CondDur;
        %CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,2);
        Factors.Position = [data(:,4), data(:,5), zeros(size(data,1),1)];
        Factors.Color = data(:,3);
    case {'CatRFfast10x10'}
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,6);
        CondTest.CondRepeat = data(:,7);
        CondTest.CondIndex = data(:,2);
        Factors.Position = [data(:,4), data(:,5), zeros(size(data,1),1)];
        Factors.Color = data(:,3);
        Factors.SpatialFreq = data(:,8);
        Factors.TemporalFreq = data(:,9);
        Factors.Ori = data(:,10);
    case {'CatRFfast'} % 'CatRFfast' is broken, condition numbers are wrong
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,6);
        CondTest.CondRepeat = data(:,7);
        % CondTest.CondIndex = data(:,2);
        Factors.Position = [data(:,4), data(:,5), zeros(size(data,1),1)];
        Factors.Color = data(:,3);
        Factors.SpatialFreq = data(:,8);
        Factors.TemporalFreq = data(:,9);
        Factors.Ori = data(:,10);
    case {'NaturalImages', 'NaturalVideos'}
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,4);
        CondTest.CondRepeat = data(:,6);
        CondTest.CondIndex = data(:,2);
        Factors.Image = data(:,2);
    case {'WholeScreenMapLP'}
        CondTest.StimOn = data(:,1);
        CondTest.StimOff = data(:,1) + data(:,6);
        CondTest.CondRepeat = data(:,7);
        CondTest.CondIndex = data(:,2);
        Factors.Position = [data(:,4), data(:,5), zeros(size(data,1),1)];
        Factors.Color = data(:,3);
    case {'CenterSurround'} % ori is -1 when either grating is off...
        CondTest.StimOn = data(:,2);
        CondTest.StimOff = data(:,2) + data(:,8);
        CondTest.CondRepeat = data(:,9);
        CondTest.CondIndex = data(:,1);
        Factors.CenterOri = data(:,6);
        Factors.Center = data(:,12) ~= -2;
        Factors.SurroundOri = data(:,12);
        Factors.Surround = data(:,12) ~= -1;
        Factors.SurroundOri(~Factors.Center) = -1;
    otherwise

        % StimTimes
        if size(data,2) < 8
            warning(['Unsupported log file: ', ex.source]);
            ex.CondTest = [];
            ex.CondRepeat = 0;
            return;
        end
        
        CondTest.StimOn = data(:,2);
        CondTest.StimOff = data(:,2) + data(:,8);
        
        Factors.SpatialFreq = data(:,3);
        Factors.TemporalFreq = data(:,4);
        Factors.Diameter = data(:,5);
        Factors.Ori = data(:,6);
        Factors.Contrast = data(:,7);
        if size(data,2) >= 10
            Factors.Velocity = data(:,10);
        else
            Factors.Velocity = data(:,4)./data(:,3);
        end
        
        % Trial no
        if size(data,2) >= 9
            CondTest.CondRepeat = data(:,9);
        end
end

% Trim CondTest
ctnames = fieldnames(CondTest);
nct = min(cellfun(@(x)length(CondTest.(x)),ctnames));
for i=1:length(ctnames)
    CondTest.(ctnames{i}) = CondTest.(ctnames{i})(1:nct);
end

% Split factors to env params and conditions
CondTestCond = [];
factorNames = fieldnames(Factors);
for f=1:length(factorNames)
    ufv = unique(Factors.(factorNames{f}));
    if length(ufv) == 1
        ex.EnvParam.(factorNames{f}) = ufv;
    else
        CondTestCond.(factorNames{f}) = Factors.(factorNames{f});
    end
end

Cond = [];
if ~isempty(CondTestCond)
    % Generate new condition numbers - old ones are sometimes wrong
    [Conditions, ~, conditionNo] = unique(struct2table(CondTestCond));
    CondTest.CondIndex = conditionNo;
    Cond = table2struct(Conditions, 'ToScalar', true);

    % Convert Cond fields to cell arrays
    CondTestCond = structfun(@(x)num2cell(x,2), CondTestCond, 'UniformOutput', false);
    Cond = structfun(@(x)num2cell(x,2), Cond, 'UniformOutput', false);
else
    CondTest.CondIndex = ones(size(data,1),1);
    CondTest.CondRepeat = 1:size(data,1);
    CondTestCond.Null = cell(1,size(data,1));
    CondTestCond.Null(:) = {0};
end

% Generate condition repeat numbers if missing
if ~isfield(CondTest, 'CondRepeat')
    repeats = zeros(max(CondTest.CondIndex),1);
    for i = 1:length(CondTest.CondIndex)
        condIndex = CondTest.CondIndex(i);
        CondTest.CondRepeat(i) = repeats(condIndex) + 1;
        repeats(CondTest.CondIndex(i)) = repeats(CondTest.CondIndex(i)) + 1;
    end
    ex.CondRepeat = max(repeats);
else
    ex.CondRepeat = max(CondTest.CondRepeat);
end

% Convert CondIndex to integers
CondTest.CondIndex = int32(CondTest.CondIndex);
CondTest.CondRepeat = int32(CondTest.CondRepeat);

ex.Cond = Cond;
ex.CondTest = CondTest;
ex.CondTestCond = CondTestCond;

end
