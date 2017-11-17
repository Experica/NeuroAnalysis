function [t] = findstatetime(stateevents,state)
%FINDSTATETIME Extract timestamp of state event
%   Detailed explanation goes here

t=[];
for i=1:length(stateevents)
    if isfield(stateevents{i},state)
        t = [t,stateevents{i}.(state)];
    end
end
if isempty(t)
    t=NaN;
end

end

