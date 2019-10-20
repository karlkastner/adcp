% Fri Feb  7 10:26:12 WIB 2014
% Karl Kastner, Berlin
%
%% converts velocity from ship to earth coordinate reference
%% expects input arguments informat:
%% vel    : float arbitrary unit
%% btvel  : float same unit as vel
%% heading_rad: float [radians]
%%
function [vel, btvel] = ship_to_earth_STATIC(vel, btvel, heading_rad, iflag)
	if (nargin() < 4)
		iflag = 0;
	end

if (0)
	% TODO verify
	% this is necessary to make the hadcp data in rasau consistent,
	% check with test_direction
	% also verified with vadcp at sanggau with datasat from 12/2013
	% note: seemed to be a regression bug, at least for the VADCP 12/2013 and the bifurcations
	% the inversion is wrong
	vel(:,:,1) = -vel(:,:,1);
	vel(:,:,2) = -vel(:,:,2);
	if (~isempty(btvel))
		btvel(:,1:2) = -btvel(:,1:2);
	end
end

	SH = sin(heading_rad);
	CH = cos(heading_rad);
	% rotate from ship reference to earth referece
	for idx=1:length(heading_rad)
	        R  = [ CH(idx) SH(idx);
                      -SH(idx) CH(idx)];
		% ship to earth is with R(H), earth to ship is with R(H)^-1,
		% see RDI coordinate transformation manual
		if (iflag)
			% invert R for ship to instrument
			% Note : R is a rotation matrix so R^-1 = R'
			R = R';
		end
	        %UV = (R*[vel(:,idx,1) vel(:,idx,2)]')';
	        UV = [vel(:,idx,1) vel(:,idx,2)]*R.';
	        vel(:,idx,1)   = UV(:,1);
	        vel(:,idx,2)   = UV(:,2);
		if (~isempty(btvel))
			btvel(idx,1:2) = btvel(idx,1:2)*R.';
			%btvel(idx,1:2) = (R*btvel(idx,1:2)')';
		end
	end % for idx
end % ship_to_earth

