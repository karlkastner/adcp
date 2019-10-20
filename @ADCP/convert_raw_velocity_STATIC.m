% Sun Jul 13 12:18:15 WIB 2014
% Karl Kastner, Berlin
%% convert raw velocity to SI units (m/s)
function [vel, btvel] = convert_raw_velocity_STATIC(vel,btvel)
	% find invalid values
	vdx = (vel   == ADCP.ERRVAL_VEL);
	bdx = (btvel == ADCP.ERRVAL_VEL);
	% transform raw value to m/s
	vel   = ADCP.VEL_SLOPE_M_PER_SEC*single(vel);
	btvel = ADCP.VEL_SLOPE_M_PER_SEC*single(btvel);
	% NaN invalid values
	vel(vdx)   = NaN;
	btvel(bdx) = NaN;
end % function convert_raw_velocity

