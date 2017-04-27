function [ p ] = lastpow2( n )
%LASTPOW2 Last higher power of 2.
%   LASTPOW2(N) returns the first P such that 2.^P <= abs(N).
%
%   Class support for input N:
%      float: double, single
%      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also LOG2, POW2, NEXTPOW2.

if ~isinteger(n)
    [f,p] = log2(abs(n));
    p = p-1;
    
    % Check for infinities and NaNs
    k = ~isfinite(f);
    p(k) = f(k);
    
else % integer case
    
end

