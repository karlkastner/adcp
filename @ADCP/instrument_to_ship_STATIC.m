% Sun Jul 13 11:21:14 WIB 2014
% Karl Kastner, Berlin
% 
%% transform velocities from instrument reference to ship reference
%% by correcting for pitch_rad and roll_rad
%%
%% input
%% vel   : float [arbitrary unit] instrument reference
%% btvel : float [arbitrary unit] instrument reference
%% pitch_rad : float [radians] true pitch_rad, not measured pitch_rad
%% roll_rad  : float [radians]
%% 
%% output
%% vel and btvel [input unit] ship reference
% 
function [vel, btvel] = instrument_to_ship_STATIC(vel,btvel,pitch_rad,roll_rad, iflag)

	if (nargin() < 5)
		iflag = 0;
	end

	CP = cos(pitch_rad);
	SP = sin(pitch_rad);
	CR = cos(roll_rad);
	SR = sin(roll_rad);

	% for each ensemble
	for idx=1:size(vel,2)
		% rotation matrix
		R =  [ 1       0        0;
	               0 CP(idx) -SP(idx);
	               0 SP(idx)  CP(idx)] ...
	           * [  CR(idx)  0 SR(idx);
	                      0  1       0;
	               -SR(idx)  0 CR(idx)];

		if (iflag)
			% invert R for ship to instrument
			% Note : R is a rotation matrix so R^-1 = R'
			R = R';
		end

		% correct all depth cells of current ping
		%vel(:,idx,1:3) = (R*squeeze(vel(:,idx,1:3))')';
		vel(:,idx,1:3) = squeeze(vel(:,idx,1:3))*R';

		% correct bottom velocity
		if (~isempty(btvel))
			btvel(idx,1:3) = btvel(idx,1:3)*R';
			%btvel(idx,1:3) = (R*btvel(idx,1:3)')';
		end
	end % for idx
end % intrument_to_ship

