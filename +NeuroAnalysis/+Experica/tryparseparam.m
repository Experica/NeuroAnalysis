function [v] = tryparseparam(name,v)
%TRYPARSEPARAM Try to convert param value to correct type
%   Detailed explanation goes here

if isa(v,'char')
    if (startsWith(v,'(') && endsWith(v,')')) || (startsWith(v,'[') && endsWith(v,']'))
        try
            v=str2double(strsplit(v(2:end-1),','));
        catch
        end
    elseif contains(v,' ')
        try
            v=str2double(strsplit(v,' '));
        catch
        end
    end
end

end