function [ result ] = EvalFun( fun,argin )
%EVALFUN Summary of this function goes here
%   Detailed explanation goes here

try
    args = NeuroAnalysis.Base.Varargin2Literal('argin',length(argin));
    eval(['result=',fun,'(',args,');']);
catch ME
    fprintf(2, '%s\n', getReport(ME));
    result.status = false;
    result.fun = fun;
    result.args = argin;
end

end

