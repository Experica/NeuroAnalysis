function [param] = standardizeParams(param, standard)
%STANDARDIZEPARAMS Summary of this function goes here
%   Detailed explanation goes here

if isempty(param)
    return;
end

import NeuroAnalysis.Base.*

% Iterate through params and make any necessary conversions
keys = fieldnames(standard);
vars = fieldnames(param);
for p = 1:length(vars) % each param
    found = 0;
    k = 1;
    while k <= length(keys) && ~found % each standard key
        type = standard.(keys{k}){1};
        if strcmp(vars{p}, keys{k}) % already standard name
            param.(keys{k}) = standardizeType(param.(keys{k}), type);
            keys = keys([1:k-1,k+1:end]);
            break;
        end
        
        options = standard.(keys{k}){2};
        o = 1;
        while o <= length(options) && ~found % each alternative
            
            if matchparam(vars{p}, options{o}, true) % case insensitive
                found = 1;
                v = getparam(param,vars{p});
                param = rmfield(param, vars{p}); % remove old name
                param.(keys{k}) = standardizeType(v, type);
                keys = keys([1:k-1,k+1:end]);
            end
            o = o + 1;
        end
        k = k + 1;
    end
end

end

function [v] = standardizeType(v, type)
% STANDARDIZETYPE convert to type if possible
    
if iscell(v)
    v = cellfun(@(x)standardizeType(x,type), v, 'UniformOutput', false);
end

switch type
    case 'double'
        if ischar(v)
            d = sscanf(v, '%f');
            if length(d) == 1
                v = d;
            end
        elseif isnumeric(v)
            if length(v) == 1
                v = double(v);
            end
        end
        
    case 'string'
        if ~ischar(v)
            v = char(string(v));
        end
        
    case 'color'
        if isnumeric(v)
            if length(v) == 1 && v > 1
                v = [repmat(double(v), 1, 3)/255, 1];
            elseif length(v) == 1
                v = [repmat(double(v), 1, 3)/255, 1];
            elseif length(v) == 3 && max(v) > 1
                v = [reshape(double(v)/255, 1, 3), 1];
            elseif length(v) == 3
                v = [reshape(double(v)/255, 1, 3), 1];
            end
        elseif ischar(v)
            switch v
                case {'Gray', 'Grey', 'gray', 'grey'}
                    v = [0.5, 0.5, 0.5, 1];
                case {'Black', 'black'}
                    v = [0, 0, 0, 1];
                case {'White', 'white'}
                    v = [1, 1, 1, 1];
            end
        end
            
    case 'bool'
        if isnumeric(v) && length(v) == 1
            v = logical(v);
        elseif ischar(v) && strcmpi(v, 'true')
            v = true;
        elseif ischar(v) && strcmpi(v, 'false')
            v = false;
        end
        
    case 'xy'
        if isnumeric(v)
            if length(v) == 2
                v = reshape(double(v), 1, 2);
            end
        elseif ischar(v)
            d1 = sscanf(v,'%f');
            d2 = sscanf(v,'%f x %f');
            if length(d1) == 2
                v = reshape(d1, 1, 2);
            elseif length(d2) == 2
                v = reshape(d2, 1, 2);
            end
        end
        
    otherwise
        warning('Unknown type %s', type);
end

end