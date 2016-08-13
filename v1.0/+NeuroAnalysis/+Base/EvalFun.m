function [ result ] = EvalFun( fun,argin )
%EVALFUN Summary of this function goes here
%   Detailed explanation goes here

eval(['result=',fun,'(',NeuroAnalysis.Base.Varargin2Literal('argin',length(argin)),');']);

end

