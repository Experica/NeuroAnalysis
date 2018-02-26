function [v] = tryparseparam(name,v)
%TRYPARSEPARAM Try to convert param value to correct type
%   Detailed explanation goes here

if isa(v,'char')
    if (startsWith(v,'(') && endsWith(v,')')) || (startsWith(v,'[') && endsWith(v,']'))
        v=str2double(strsplit(v(2:end-1),','));
    end
end

end

