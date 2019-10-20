% Sun Jul 13 12:30:09 WIB 2014
% Karl Kastner, Berlin
%
%% convert raw instrument heading angle to [rad]
function heading_rad = heading_rad(obj)
	% raw values to radians
	heading_rad = deg2rad(ADCP.HEADING_SLOPE * single(obj.dat.heading));
	% validate heading
	% TODO also validate roll and pitch (note, roll and pitch are signed)
	fdx = obj.dat.heading > ADCP.HEADING_MAX;
	if (sum(fdx) > 0)
		disp('waring: Implausible heading provided by instrument');
		heading_rad(fdx) = NaN;
	end
end

