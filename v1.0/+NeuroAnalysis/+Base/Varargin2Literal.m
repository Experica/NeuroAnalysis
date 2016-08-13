function [ vail ] = Varargin2Literal( vainame,vain )
%VARARGIN2LITERAL Summary of this function goes here
%   Detailed explanation goes here

v=arrayfun(@(i)[vainame,'{',num2str(i),'}'],1:vain,'UniformOutput',false);
vail = strjoin(v,',');
end

