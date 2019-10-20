% Sun Jul 13 11:08:40 WIB 2014
% Karl Kastner, Berlin

%% transform the 4 beam velocities into a set of 3 orthogonal velocities
%% and 1 error velocity
%%
%% input
%% vel       : float [arbitrary unit] beam reference system
%% btvel     : float [arbitrary unit] beam reference system
%% beamangle : float [radians]
%%
%% output
%% vel and btvel [input unit] instrument reference system
%%
%% TODO account for NaNs either by three beam solution or interpolation
%%
function [vel, btvel] = beam_to_instrument_STATIC(vel,btvel,beamangle_rad, mode, iflag)
	if (nargin() < 3 || isempty(beamangle_rad))
		beamangle_rad = VADCP.DEFAULT_BEAMANGLE_RAD;
	end

	% convex transducer head (c = +1)
	% TODO, make this and argument
	c = +1;
	a = 1/(2*sin(beamangle_rad));
	b = 1/(4*cos(beamangle_rad));
	d = a/sqrt(2);
	% transformation matrix (same for all ensembles)
	T = [ c*a -c*a    0   0;
		0    0 -c*a c*a;
		b    b    b   b;
		d    d   -d  -d ];
	if (nargin() > 4 && ~isempty(iflag) && iflag)
		T = pinv(T);
	end
	rows = size(vel,1);
	cols = size(vel,2);
	vel  = reshape(vel,rows*cols,4);
	%vel  = (T*vel')';
	vel  = vel*T';
	vel = reshape(vel,[rows cols 4]);
	if (~isempty(btvel))
		%btvel = (T*btvel')';
		btvel = btvel*T';
	end
end % beam_to_instrument_STATIC

