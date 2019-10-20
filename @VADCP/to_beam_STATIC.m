% Thu Jul 17 18:59:33 WIB 2014
% Karl Kastner, Berlin
%% convert velocity data to beam reference
function [vel, btvel] = to_beam_STATIC(corstr,FileNumber,vel,btvel,heading,pitch,roll)
	for idx=1:length(corstr)
		fdx = find(FileNumber == idx);
		inverse = true;
		flag = false;
		if (~isempty(btvel))
			if (strcmp(corstr{idx}(1:4),'Earth'))
				[vel.ship(:,fdx,:), btvel.ship(fdx,:)] = ADCP.ship_to_earth_STATIC(vel.earth(:,fdx,:),btvel.earth(fdx,:),heading(fdx),inverse);
				flag = true;
			end
			if (strcmp(corstr{idx}(1:4),'Ship') || flag)
				% TODO roll minus pi ? is lacking in the SHIP case
				[vel.inst(:,fdx,:), btvel.inst(fdx,:)] = ADCP.instrument_to_ship_STATIC(vel.ship(:,fdx,:),btvel.ship(fdx,:),pitch(fdx),roll(fdx),inverse);
				flag = true;
			end
			if (strcmp(corstr{idx}(1:4),'Inst') || flag)
				[vel.beam(:,fdx,:), btvel.beam(fdx,:)] = VADCP.beam_to_instrument_STATIC(vel.inst(:,fdx,:),btvel.inst(fdx,:),[],[],inverse);
				flag = true;
			end
			if (strcmp(corstr{idx}(1:4),'Beam') || flag)
				% nothing to do
			else
				error(['Unknown coordinate system ' corstr{idx}]);
			end
		else
			if (strcmp(corstr{idx}(1:4),'Earth'))
				[vel.ship(:,fdx,:)] = ADCP.ship_to_earth_STATIC(vel.earth(:,fdx,:),[],heading(fdx),inverse);
				flag = true;
			end
			if (strcmp(corstr{idx}(1:4),'Ship') || flag)
				% TODO roll minus pi ? is lacking in the SHIP case
				[vel.inst(:,fdx,:)] = ADCP.instrument_to_ship_STATIC(vel.ship(:,fdx,:),[],pitch(fdx),roll(fdx),inverse);
				flag = true;
			end
			if (strcmp(corstr{idx}(1:4),'Inst') || flag)
				[vel.beam(:,fdx,:)] = VADCP.beam_to_instrument_STATIC(vel.inst(:,fdx,:),[],[],[],inverse);
				flag = true;
			end
			if (strcmp(corstr{idx}(1:4),'Beam') || flag)
				% nothing to do
			else
				error(['Unknown coordinate system ' corstr{idx}]);
			end
		end % isempty btvel
	end % for idx
end % to_earth()

