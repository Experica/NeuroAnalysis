function [ex] = StandardizeEx(ex)
%STANDARDIZEEX Summary of this function goes here
%   Detailed explanation goes here

import NeuroAnalysis.Base.*

%% Experimental data
ex = standardizeParams(ex, ExStandard);

%% Environment parameters and factors
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

