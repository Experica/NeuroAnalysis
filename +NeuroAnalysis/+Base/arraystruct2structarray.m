function [sa] = arraystruct2structarray(as)
%ARRAYSTRUCT2STRUCTARRAY Summary of this function goes here
%   Detailed explanation goes here

fs = fieldnames(as);
sa=struct;
for i = 1:length(fs)
    sa.(fs{i})=[as.(fs{i})];
end

end

