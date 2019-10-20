% Sun Jul 13 12:30:09 WIB 2014
% Karl Kastner, Berlin
%
%% convert raw instrument roll angle to [rad]
%
function roll_rad = roll_rad(obj)
	roll_rad = deg2rad(ADCP.ROLL_SLOPE  * single(obj.dat.roll));
	% adapt roll for downward looking deployment
	if (nargin() < 5 || 0 == upward)
		roll_rad = pi - roll_rad;
	end
end

