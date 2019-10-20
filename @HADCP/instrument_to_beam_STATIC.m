% Mon Jul 14 01:43:06 WIB 2014
% Karl Kastner, Berlin

%% transform the 3 beam velocities into a set of 2 orthogonal velocities
%% and 1 error velocity
%% This uses always three beams (no two beam solutions)
%%
%% input
%% vel       : float [arbitrary unit] beam reference system
%% btvel     : float [arbitrary unit] beam reference system
%% beamangle : float [radians]
%%
%% output
%% vel and btvel [input unit] instrument reference system
%%
%% mode : beams used for all transformations
%%        123, 12, 23, 13
% TODO account for NaNs either by three beam solution or interpolation
%
function [vel btvel] = beam_to_instrument(vel,btvel,beamangle_rad,mode,iflag)
	if (nargin() < 4)
		mode = 123;
	end
	if (nargin() < 5)
		iflag = 0;
	end

	% convex transducer head (c = +1)
	a  = 1/(2*sin(beamangle_rad));
	b  = cos(beamangle_rad)/(1+2*cos(beamangle_rad)^2);
	c  = 1/(1+2*cos(beamangle_rad)^2);
	e1 = 1/(2*sin(beamangle_rad)*sqrt(1+2*cos(beamangle_rad)^2));
	e2 = 1/(  tan(beamangle_rad)*sqrt(1+2*cos(beamangle_rad)^2));
	% transformation matrix (same for all ensembles)
	% note : column and row 4 are all zeros, 
	% so left out, this also yields the only valid result if there are
	% nans in the redundant beam
	switch (mode)
		case {123}
			T = [   -a     a     0  0;
				-b    -b    -c  0;
                                 0     0     0  0
				e1    e1   -e2  0];
		case {12}
			T = [  -a      a      0 0;
                               -e1/e2  -e1/e2 0 0;
                                0      0      0 0;
                                0      0      0 0];
		case {13}
			% TODO, this was not given, made up by myself
			T = [ -2*a 0  a*e2/e1 0;
                               0   0       -1 0;
                               0   0        0 0;
                               0   0        0 0];
		case {23}     
			T = [ 0 2*a -a*e2/e1 0;
                              0   0       -1 0;
                              0   0        0 0;
                              0   0        0 0];
		otherwise
			error('unknown beam combination');
	end
	if (iflag)
		T = pinv(T);
	end
	rows = size(vel,1);
	cols = size(vel,2);
	vel  = reshape(vel,rows*cols,4);
	vel  = (T*vel(:,1:4)')';
	vel = reshape(vel,[rows cols 4]);
%	vel(:,1:3)  = (T*vel(:,1:3)')';
	%vel(:,:,1:3) = reshape(vel(:,1:3),[rows cols 3]);
%	vel = reshape(vel(:,1:3),[rows cols 3]);
%	vel(:,:,4) = 0;
	if (~isempty(btvel))
		btvel(:,1:4) = (T*btvel(:,1:4)')';
	end
	%	btvel(:,4) = 0;
end % beam_to_instrument

