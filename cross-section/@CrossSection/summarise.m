% Mon 30 May 13:59:54 CEST 2016
% Karl Kastner, Berlin
%% summarize discharge of cross section
function [discharge, area] = summarize(obj)
	discharge_A     = obj.discharge;
	area_A          = obj.area;
	discharge.total = mean(discharge_A.total);
	area.total      = mean(area_A.total);
end

