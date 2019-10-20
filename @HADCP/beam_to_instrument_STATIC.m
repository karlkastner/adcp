% Mon Jul 14 01:43:06 WIB 2014
% Karl Kastner, Berlin
%
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
%
% TODO account for NaNs either by three beam solution or interpolation
%
% mode : beams used for all transformations
%        123, 12, 23, 13
function [ovel btvel] = beam_to_instrument_STATIC(vel,btvel,beamangle_rad,mode,iflag)
	if (nargin() < 4 || isempty(mode))
		mode = 'auto';
	end
	if (nargin() < 5)
		iflag = 0;
	end

	nrow = size(vel,1);
	ncol = size(vel,2);
	n    = nrow*ncol;
	vel  = reshape(vel,n,4);
	ovel = NaN(size(vel));

	% swap, for overhead deployment
	swapflag = true;
	if (swapflag) % && ~iflag)
%		vel(:,1:2) = vel(:,1:2)*[0 1; 1 0];
%		vel(:,1:2) = -vel(:,1:2);
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
	T123 = [   -a     a     0  0;
   		   -b    -b    -c  0;
                    0     0     0  0;
		    e1    e1   -e2 0];
	T12 = [  -a      a      0 0;
               -e1/e2  -e1/e2   0 0;
                  0      0      0 0;
                  0      0      0 0];
	% These are not given in the manual, and derived manually
	T13 = [ -2*a 0  a*e2/e1 0;
                 0   0       -1 0;
                 0   0        0 0;
                 0   0        0 0];
	T23 = [ 0 2*a -a*e2/e1 0;
                0   0       -1 0;
                0   0        0 0;
                0   0        0 0];

	switch (mode)
		case {12}
			f = 3*ones(n,1);
		case {13}
			f = 5*ones(n,1);
		case {23}
			f = 6*ones(n,1);
		case {123}
			f = 7*ones(n,1);
		case {'auto'}
			if (iflag)
				% This is not unique, as the ADCP does not store which beams are used
				% assume 3-beam solution
				f = 7*ones(n,1);
			else
				f = isfinite(vel(:,1)) + 2*isfinite(vel(:,2)) + 4*isfinite(vel(:,3));
			end
		otherwise
			error('unknown beam combination');
	end
	if (iflag)
		T12  = pinv(T12);
		T13  = pinv(T13);
		T23  = pinv(T23);
		T123 = pinv(T123);
	end
	% NaN to zeros, otherwise as 0*NaN is NaN
	fdx  = isnan(vel);
	vel(fdx) = 0;

	% 12
	fdx = (f==3);
	ovel(fdx,:) = vel(fdx,:)*T12';
	% 13
	fdx = (f==5);
	ovel(fdx,:) = vel(fdx,:)*T13';
	% 23
	fdx = (f==6);
	ovel(fdx,:) = vel(fdx,:)*T23';
	% 123
	fdx = (f==7);
	ovel(fdx,:) = vel(fdx,:)*T123';

	ovel = reshape(ovel,[nrow ncol 4]);
	if (~isempty(btvel))
		btvel = btvel*T123';
	end

%	if (swapflag && iflag)
%		vel(:,1:2) = vel(:,1:2)*[0 1; 1 0];
%		vel(:,1:2) = vel(:,1:2)
%	end
end % beam_to_instrument

