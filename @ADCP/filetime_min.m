% Thu 20 Jul 17:10:39 CEST 2017
%% start time of each file
function [t, id, obj] = filetime_min(obj)
	fdx = find([1; cvec(diff(obj.dat.FileNumber))]);
	id = obj.dat.FileNumber(fdx);
	t = obj.time(fdx);
end
