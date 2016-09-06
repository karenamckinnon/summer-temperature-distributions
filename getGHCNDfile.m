function [T,d,flag] = getGHCNDfile(yr_begin,yr_end,varname,stationRecord)
% [T,dates,flag] = get_data(yr_begin,yr_end,varname,stationRecord)
% Extracts data for a given variable from the stationRecord structure, 
% subsetting for the desired time period.
%     
% In:  
%    yr_begin - 
%    yr_end   - 
%    varname  - 
%    stationRecord - 
% Out:
%    T - the data for varname
%    dates - the dates for those data, in matlab datenum format
%    flag  - returns zero if the record doesn't contain data
%            matching the input criteria.

flag = 1; % set to zero if record doesn't fit date criteria
d = stationRecord.dates; 
T = 0;
% only consider stations that have record until specified yr_end ...
if max(year(d)) >= yr_end
    % ... and that begin earlier than yr_begin
    % set yr_begin = 0 if no preference on earlier year
    if min(year(d)) <= yr_begin
        yr_range = yr_begin:yr_end;

        try % Some weird files have numel(yr_range) = 0
            t1 = find(year(d) == yr_range(1),1,'first');
            t2 = find(year(d) == yr_range(end),1,'last');

            fn = char(['stationRecord.' varname '.data']);
            T = eval(fn);
            T = double(T(t1:t2));
            d = d(t1:t2);

        catch err
            flag = 0;
        end

    else
        flag = 0;
    end
else
    flag = 0;
end

return