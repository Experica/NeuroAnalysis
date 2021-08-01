function [ex] = StandardizeEx(ex)
%STANDARDIZEEX Convert params to standard name and value type.
%   Detailed explanation goes here

import NeuroAnalysis.Base.*

%% Experimental parameters
ex = standardizeParams(ex, ExStandard);

%% Environmental parameters and factors
if isfield(ex, 'EnvParam')
    ex.EnvParam = standardizeParams(ex.EnvParam, ExFactorStandard);
end
if isfield(ex, 'Cond')
    ex.Cond = standardizeParams(ex.Cond, ExFactorStandard);
end
if isfield(ex, 'CondTestCond')
    ex.CondTestCond = standardizeParams(ex.CondTestCond, ExFactorStandard);
end

end

