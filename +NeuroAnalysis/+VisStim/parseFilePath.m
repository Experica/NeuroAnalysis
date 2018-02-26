% --- parse a filename
function [animalID, fileNo, unit, stimType] = parseFilePath (filePath)
% FileName  should be a string of the format AnimalID#FileNo[StimName]
%
% AnimalID  string for the animal ID
% FileNo    integer for the file number
% Unit      string for recording site / unit
% StimName  string describing the stimulus choice

animalID = '';
fileNo = NaN;
stimType = '';

pathParts = strsplit(filePath, filesep);
unit = char(pathParts{end-1});

fileName = strrep(pathParts{end}, ']', ' ');
[nums, c] = sscanf(fileName, '%c%4d#%d[%s*');
if c < 4; return; end
animalID = [char(nums(1)), num2str(nums(2))];
fileNo = nums(3);
stimType = char(nums(4:end))';