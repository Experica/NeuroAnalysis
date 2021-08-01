classdef ExFactorStandard
    %EXDATASTANDARD Case-insensitive naming conventions
    %   Detailed explanation goes here
    
    properties
        
        Ori = {'double', {'ori', 'orientation'}}
        Diameter = {'double', {'diameter', 'diam', 'aperture', 'apt'}}
        Size = {'xy', {'size'}}
        Color = {'color', {'color'}}
        ColorID = {'int', {'colormod'}}
        Visible = {'bool', {'visible', 'on'}}
        Position = {'xy', {'position', 'xyposition', 'location'}}
        Luminance = {'double',{'luminance'}}
        Contrast = {'double', {'contrast', 'con'}}
        SpatialFreq = {'double', {'spatialfreq', 'spatialfrequency', 'sf', 'sfreq','s_freq'}}
        TemporalFreq = {'double', {'temporalfreq', 'temporalfrequency', 'tf', 'tfreq'}}
        SpatialPhase = {'double',{'spatialphase'}}
        ModulateTemporalFreq = {'double',{'modulatetemporalfreq'}}
        GratingType = {'string', {'gratingtype', 'grating'}}
        BGColor = {'color', {'bgcolor', 'bg', 'blankcolor'}}
        MinColor = {'color', {'mincolor'}}
        MaxColor = {'color', {'maxcolor'}}
        MaskRadius = {'double',{'maskradius'}}
        MaskSigma = {'double',{'masksigma'}}
        MaskType = {'string',{'masktype'}}
        ModulateType = {'string',{'modulatetype'}}
        Duty = {'double',{'duty'}}
        ModulateDuty = {'double',{'modulateduty'}}
        
    end
    
end

