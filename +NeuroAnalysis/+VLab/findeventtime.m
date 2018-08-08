function [t] = findeventtime(eventtimes,event)
%FINDEVENTTIME Find timestamps of event
%   Detailed explanation goes here

t=[];
for i=1:length(eventtimes)
    if isfield(eventtimes{i},event)
        t=[t,eventtimes{i}.(event)];
    end
end

end

