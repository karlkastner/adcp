% Tue Jan 28 09:29:22 WIB 2014
% Karl Kastner, Berlin
%
%% fill gaps in ensemble coordinates
%
% TODO introduce a separate fields Xi, Yi for fixed coordinates or compute this on the fly
% TODO use btrack to fill gaps
% TODO warn if not sorted
% TODO use vtg
% TODO combine gps and bottom track with a kalman filter
function obj = fill_coordinate_gaps(obj)
	time = obj.time;
	X    = obj.X;
	Y    = obj.Y;

	% sort samples
	time0 = time;

	[time, sdx] = sort(time);
	X = X(sdx);
	Y = Y(sdx);

	dtmax     = obj.dt_max_gap; 
	X         = fixnan(time,X,dtmax);
	Y         = fixnan(time,Y,dtmax);

	obj.X(sdx) = X;
	obj.Y(sdx) = Y;
end % fill_coordinate_gaps

