% Sun Jul 13 20:14:12 WIB 2014
% Karl Kastner, Berlin
%% wrappter to transform velocities to world coordinate reference
function [vel, btvel] = to_earth_STATIC(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad,beamangle_rad,mode)
	if (nargin() < 8)
		beamangle_rad = HADCP.DEFAULT_BEAMANGLE_RAD;
	end
	if (nargin() < 9)
		mode = [];
	end
	for idx=1:length(corstr)
		fdx = find(FileNumber == idx);
		if (~isempty(btvel))
			btvel_i = btvel(fdx,:);
		else
			btvel_i = [];
		end
		flag = 0;
		if (strcmp(corstr{idx}(1:4),'Beam'))
				[vel.inst(:,fdx,:) btvel_i] = HADCP.beam_to_instrument_STATIC(vel.beam(:,fdx,:),btvel_i,beamangle_rad(idx),mode);
				flag = 1;
		end
		if (flag | strcmp(corstr{idx}(1:4),'Inst'))
				[vel.ship(:,fdx,:) btvel_i] = ADCP.instrument_to_ship_STATIC(vel.inst(:,fdx,:),btvel_i,pitch_rad(fdx),roll_rad(fdx));
				flag = 1;
		end
		if (flag | strcmp(corstr{idx}(1:4),'Ship'))
				[vel.earth(:,fdx,:) btvel_i] = ADCP.ship_to_earth_STATIC(vel.ship(:,fdx,:),btvel_i,heading_rad(fdx));
		end
		if (~isempty(btvel))
			btvel(fdx,:) = btvel_i;
		end
	end % for idx

end % to_earth()

