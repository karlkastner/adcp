% Sun Jul 13 12:30:09 WIB 2014
% Karl Kastner, Berlin
%
%% convert raw pitch to radians
%
function pitch_rad = pitch_rad(obj)
		% instrument measured pitch angle [rad] (not true pitch)
		pitch_rad = deg2rad(ADCP.PITCH_SLOPE * single(obj.dat.pitch));
		roll_rad = obj.roll_rad();

		% adapt pitch as pitch and roll not measured by gimbals (affecting each other)
		pitch_rad = atan(tan(pitch_rad).*cos(roll_rad));
end % pitch_rad

