classdef ExFactorStandard
    %EXDATASTANDARD Case-insensitive naming conventions
    %   Detailed explanation goes here
    
    properties
        
        Ori = {'double', {'ori', 'orientation'}}
        Diameter = {'double', {'diameter', 'diam', 'aperture', 'apt'}}
        Size = {'xy', {'size'}}
        Color = {'color', {'color'}}
        Visible = {'bool', {'visible', 'on'}}
        Position = {'xy', {'position', 'xyposition', 'location'}}
        Contrast = {'double', {'contrast', 'con'}}
        SpatialFreq = {'double', {'spatialfreq', 'spatialfrequency', 'sf', 'sfreq','s_freq'}}
        TemporalFreq = {'double', {'temporalfreq', 'temporalfrequency', 'tf', 'tfreq'}}
        GratingType = {'string', {'gratingtype', 'grating'}}
        BGColor = {'color', {'bgcolor', 'bg', 'blankcolor'}}
        
    end
    
end

