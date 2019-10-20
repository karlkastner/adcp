% Do 5. Nov 14:45:11 CET 2015
% Karl Kastner, Berlin
%% split data set into specific time slots
% TODO this is crappy, time slot splitting should not depend on function or other modes
%
function [obj] = determine_time_slots(obj,adcp)
	tmax                = max(adcp.time);
	tmin 		    = min(adcp.time);

%	if (strcmp('FUNCTION',upper(obj.tmode)))
	%fdx                 = isfinite(cvec(obj.level.time).*cvec(obj.level.val));
	%adcp.ens.level      = interp1(obj.level.time(fdx),obj.level.val(fdx),adcp.time,'linear');


	% time for discharge computation
	% TODO tmode FUNCTION is deprecated
%	if (strcmp('FUNCTION',upper(obj.tmode)))
	dt = obj.dt;
	if (isempty(dt))
		n = 1;
	else
		n = max(1,round( (tmax-tmin)/dt ));
	end
	obj.t0		    = tmin + (tmax-tmin)*(0:n)/(n);
%	else
%		adcp.ens.level      = zeros(adcp.ens.n,1,'single');
%		adcp.ens.level4     = zeros(adcp.ens.n,4,'single');
		%n                   = 1;
%		obj.t0		    = tmin;
%	end
	%if (strcmp('FUNCTION',upper(obj.tmode)))
	%obj.level0          = obj.levelt(obj.t0);
	%end
end % determine_time_slots

