%
%copies a list of fields from a source to a target structure
%
%fields -> cell array of strings. entries in the cell array can also be an other cell array
%with two string entries: src and target. (to rename fields).
%
%i.e. fieldsSrc = {'a1',{'a2','a3',fun}, {'from','to',@int32})
%
%will copy target.a1=src1 and target.a3=fun(src2).
%
%if fieldsSrc is not supplied, copies all fields in src to target structure.
%
%
%urut/may07
%modified leo scholl december 2017
function target = copyStructFields(src,target,fieldsSrc,valueFun)
if ~isstruct(src)
    warning('src is not a struct - ignore. nothing is copied.');
    return;
end

if nargin==2 || isempty(fieldsSrc)
    fieldsSrc=fieldnames(src);
end

if nargin<4
    valueFun = [];
end

for i=1:length(fieldsSrc)
    fieldSrc=fieldsSrc{i};
    
    if iscell(fieldSrc)
        if length(fieldSrc) > 2
            valueFun=fieldSrc{3};
        end
        fieldTarget=fieldSrc{2};
        fieldSrc=fieldSrc{1};
    else
        fieldTarget=fieldSrc;
    end
    
    if isfield(src,fieldSrc) 
        if ~isempty(valueFun)
            target.(fieldTarget) = valueFun(src.(fieldSrc));   %a dynamic field name instead of eval
        else
            target.(fieldTarget) = src.(fieldSrc);
        end
    end
end