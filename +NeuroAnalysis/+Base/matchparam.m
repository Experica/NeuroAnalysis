function [ismatch] = matchparam(param,str,ignorecase)
%MATCHPARAM Summary of this function goes here
%   Detailed explanation goes here

if nargin==2
    ignorecase=false;
end

ismatch = false;
if ignorecase && any((strcmpi(param, str)) || ...
        any(startsWith(param, [str, '0x40'], 'IgnoreCase', true)))
    ismatch = true;
elseif any(strcmp(param, str)) || any(startsWith(param,[str, '0x40']))
    ismatch = true;
end

end

