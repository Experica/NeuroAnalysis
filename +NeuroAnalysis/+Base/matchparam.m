function [match] = matchparam(param,str,ignorecase)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
match = false;
if ignorecase && (strcmpi(param, str) || ...
        startsWith(param, [str, '0x40'], 'IgnoreCase', true))
    match = true;
elseif strcmp(param, str) || startsWith(param,[str, '0x40'])
    match = true;
end
end

