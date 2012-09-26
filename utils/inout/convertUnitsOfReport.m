function report = convertUnitsOfReport(report)
% convert report object to standard metric units
%report = convertUnitsOfReport(report)
report.TIME =report.TIME/day;
field_names = {'WWPR','WOPR','WVPT'};
for i=1:numel(field_names)
   name = char(field_names(i));
   report.(name) = report.(name)*day;
end
report.WBHP = report.WBHP/barsa;