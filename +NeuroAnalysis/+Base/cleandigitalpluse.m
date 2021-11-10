function [time,data] = cleandigitalpluse(time,data,width,sign,info)
%CLEANDIGITALPLUSE remove narrow digital pluses

if nargin < 3
    width = 1;
    sign = true;
    info = '';
elseif nargin < 4
    sign = true;
    info = '';
elseif nargin < 5
    info = '';
end

if sign % positive pluse
    ri = find(data~=0);
    fi = ri+1;
    if fi(end)>length(data)
        fi(end) = [];
        ri(end) = [];
    end
    i = find(time(fi)-time(ri) < width);
else % negative pluse
    fi = find(data==0);
    ri = fi+1;
    if ri(end)>length(data)
        fi(end) = [];
        ri(end) = [];
    end
    i = find(time(ri)-time(fi) < width);
end

if ~isempty(i)
    if sign
        warning('Clean Positive Digital Pluse Width < %g, At _|%s    %s|_ %s ...',width,num2str(ri(i)),num2str(fi(i)),info);
    else
        warning('Clean Negative Digital Pluse Width < %g, At -|%s    %s|- %s ...',width,num2str(fi(i)),num2str(ri(i)),info);
    end
    time([ri(i),fi(i)])=[];
    data([ri(i),fi(i)])=[];
end

end

