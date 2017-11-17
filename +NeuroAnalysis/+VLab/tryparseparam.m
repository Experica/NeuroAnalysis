function [v] = tryparseparam(name,v)
%TRYPARSEPARAM Try to get correct type of param value
%   Detailed explanation goes here

    function v=parsevaluestring(v)
        if (startsWith(v,'(') && endsWith(v,')')) || (startsWith(v,'[') && endsWith(v,']'))
            v=str2double(strsplit(v(2:end-1),','));
        end
    end

    function v=parseparam(v)
        if isa(v,'char')
            v=parsevaluestring(v);
        end
    end

switch name
    otherwise
        v=parseparam(v);
end

end

