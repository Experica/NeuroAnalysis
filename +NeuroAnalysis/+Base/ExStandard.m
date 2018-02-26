classdef ExStandard
    %EXDATASTANDARD Case-insensitive naming conventions
    %   Detailed explanation goes here
    
    properties
        
        Subject_ID = {'string', {'subject_id', 'subjectid', 'animalid'}}
        RecordSession = {'string', {'recordingsession', 'recordsession'}}
        RecordSite = {'string', {'recordingsite', 'recordsite', 'unit'}}
        ID = {'string', {'stimtype', 'id', 'stimid'}}
        File_ID = {'string', {'file_id', 'fileid', 'fileno'}}
        
    end

end

