% Mon Jul 14 00:02:55 WIB 2014
% Karl Kastner, Berlin
%
%% convert raw pressure to bar
%
function pressure_bar = pressure_bar(obj)
%function [pressure_bar, idepth_m] = convert_raw_pressure_STATIC(pressure_raw)
	pressure_raw = single(obj.dat.pressure);

	% find errornous values
	fdx = (pressure_raw > ADCP.PRESSURE_ERR_THRESH);

	% convert raw pressure value into bar
	% pressure had more than 24 significant bits, double is necessary
	pressure_bar = ADCP.PRESSURE_SLOPE_BAR*double(pressure_raw) + ADCP.PRESSURE_OFFSET_BAR;

	% invalidate error values
	pressure_bar(fdx) = NaN;
	
end % function pressure_bar

