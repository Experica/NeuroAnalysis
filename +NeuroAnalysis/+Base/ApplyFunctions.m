function [ results ] = ApplyFunctions( funlist,vararginlist,isparallel )
%APPLYFUNCTIONS Apply list of function to list of arguments, and return list of result.
%   Since matlab doesn't support multiple value return as tuple,
%   only the first value is captured in result list.

if nargin==2
    isparallel = false;
end

funn = length(funlist);
argn = length(vararginlist);
if funn~=argn
    error('Number of function (%d) does not match number of arguments (%d).',funn,argn);
end

if isparallel
    parfor i=1:funn
        fun = funlist{i};
        argin = vararginlist{i};
        disp(['Parallel Apply Function:    ',fun,'    ===========================================>']);
        results{i} = NeuroAnalysis.Base.EvalFun(fun,argin);
    end
else
    for i=1:funn
        fun = funlist{i};
        argin = vararginlist{i};
        disp(['Apply Function:    ',fun,'    ===========================================>']);
        results{i} = NeuroAnalysis.Base.EvalFun(fun,argin);
    end
end

end

