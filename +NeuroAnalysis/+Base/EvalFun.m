function [ result ] = EvalFun( fun,argin )
%EVALFUN Summary of this function goes here
%   Detailed explanation goes here

try
    args = NeuroAnalysis.Base.Varargin2Literal('argin',length(argin));
    eval(['result=',fun,'(',args,');']);
catch ME
    warning([ME.identifier,': ',ME.message]);
    for i=1:length(ME.stack)
        ME.stack(i);
    end
    result.status = false;
    result.fun = fun;
    result.args = args;
end

end

