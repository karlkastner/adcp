% Sun 27 May 13:09:33 CEST 2018
%% assign water level to adcp ensembles (combine gauge with boat data)
function obj = assign_water_level(obj,time,level)
	% TODO save
	level              = interp1(time,level,obj.time);
	obj.ens.level      = cvec(level);
end

