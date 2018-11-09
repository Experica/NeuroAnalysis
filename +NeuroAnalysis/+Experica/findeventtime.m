function [ts] = findeventtime(eventtimes,events)
%FINDEVENTTIME Find timestamp of events
%   Detailed explanation goes here

ts=[];
sidx=1;
net=length(eventtimes);
for j=1:length(events)
    for i=sidx:net
        if isfield(eventtimes{i},events{j})
            ts=[ts,eventtimes{i}.(events{j})];
            sidx=i+1;
            break;
        end
    end
end

end