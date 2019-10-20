% Sat Jan 10 20:38:23 CET 2015
% Karl Kastner, Berlin
%% get the cartesian (world) coordinates of the HADCP central beam bins
% TODO this should be part of the HADCP class and do a different implementation for the VADCP
% TODO quick hack for HADCP bin range - central beam(vel distances) not scaled by 1/cos
function [X,Y,obj] = calc_bin_coordinates(obj)

	% heading (assumed constant over deployment)
	% TODO allow for variable heading
	%mheading_rad = nanmedian(obj.heading_rad);
	mheading_rad = meanangle(obj.heading_rad);

	% TODO, do this for all beams
	% angle is with respect to to North, therefore dx ~ sin and dy ~ cos
	X = double(obj.R*sin(mheading_rad));
	Y = double(obj.R*cos(mheading_rad));
	if (~isempty(obj.X))
		X = obj.X + X;
		Y = obj.Y + Y;
	else
		%X = X;
		%Y = Y;
		warning('No instrument coordinates given, coordinates are with respect to origin of coordinate system');
	end
end % calc_bin_coordinates

