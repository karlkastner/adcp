% Sat Jan 25 13:35:28 WIB 2014
% Karl Kastner, Berlin
%
%% claculate positions in world coordinates where the individual beams hit the bottom
%
% TODO, this does not account for pitch and roll
function obj = calc_beamcoords(obj)

	if (~isempty(obj.adcp.btrange))

	% expand X and Y coordinate for each beam
	if (~isempty(obj.X))
		X4 = repmat(obj.X,1,4);
		Y4 = repmat(obj.Y,1,4);
	else
		X4 = 0;
		Y4 = 0;
	end

	% correct X and Y position of each beam for attitude angle
	% determine xy-position of beam intersection with the bottom
	% x and y coordinates of the beam intersection
	% beams are arranged : 3 = 0, 1 = 90, 4 = 180, 2 = 270
	% 0 is north
	beamangle_rad = obj.adcp.beamangle_rad;
	beamangle_rad = beamangle_rad(1);
	heading_rad   = cvec(obj.adcp.heading_rad);

	R    = tan(beamangle_rad)*obj.adcp.btrange;
	sinh = sin(heading_rad(:));
	cosh = cos(heading_rad(:));

	% transducers are arranged 3(0),1(90),4(180),2(270)
	% 1 and 2 swapped when downward looking
	%  downward
	%  b dir     x  y
	%  3 north   s  c
	%  2 east    c -s
	%  4 south  -s -c
	%  1 west   -c  s
	dx   = [-cosh,  cosh,  sinh, -sinh].*R;
	dy   = [ sinh, -sinh,  cosh, -cosh].*R;

	X4 = X4 + double(dx);
	Y4 = Y4 + double(dy);

	% write back
	obj.X4 = X4;
	obj.Y4 = Y4;

	end
end % calc_beamcoords

