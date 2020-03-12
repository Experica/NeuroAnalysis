function [ts] = findeventtime(eventtimes,events)
%FINDEVENTTIME Find timestamp of events
%   Detailed explanation goes here

ts=[];
start=1;
net=length(eventtimes);
for j=1:length(events)
    for i=start:net
        if isfield(eventtimes{i},events{j})
            ts=[ts,eventtimes{i}.(events{j})];
            start=i+1;
            break;
        end
    end
end

end