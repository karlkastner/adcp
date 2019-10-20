% Thu Jul 17 18:59:33 WIB 2014
% Karl Kastner, Berlin
%% transform coordinates to cartesian world reference system (earth)
function [vel, btvel] = to_earth_STATIC(corstr,FileNumber,vel,btvel,heading,pitch,roll)
	for idx=1:length(corstr)
		fdx = find(FileNumber == idx);
		if (~isempty(btvel))
		switch (corstr{idx}(1:4))
		case {'Beam'}
			[vel.inst(:,fdx,:),  btvel.inst(fdx,:)]  = VADCP.beam_to_instrument_STATIC(vel.beam(:,fdx,:),btvel.beam(fdx,:));
			[vel.ship(:,fdx,:),  btvel.ship(fdx,:)]  = ADCP.instrument_to_ship_STATIC(vel.inst(:,fdx,:),btvel.inst(fdx,:),pitch(fdx),roll(fdx)-pi);
			[vel.earth(:,fdx,:), btvel.earth(fdx,:)] = ADCP.ship_to_earth_STATIC(vel.ship(:,fdx,:),btvel.ship(fdx,:),heading(fdx));
			disp('Beam coordinates found');
		case {'Inst'}
			[vel.ship(:,fdx,:),  btvel.ship(fdx,:)]  = ADCP.instrument_to_ship_STATIC(vel.inst(:,fdx,:),btvel.inst(fdx,:),pitch(fdx),roll(fdx)-pi);
			[vel.earth(:,fdx,:), btvel.earth(fdx,:)] = ADCP.ship_to_earth_STATIC(vel.ship(:,fdx,:),btvel.ship(fdx,:),heading(fdx));
			disp('Instrument coordinates found');
		case {'Ship'}
			[vel.earth(:,fdx,:), btvel.earth(fdx,:)] = ADCP.ship_to_earth_STATIC(vel.ship(:,fdx,:),btvel.ship(fdx,:),heading(fdx));
			disp('Ship coordinates found');
		case {'Eart'}
			disp('Earth coordinates found');
			% nothing to do
		otherwise
			error(['Unknown coordinate system ' corstr{idx}]);
		end % switch

		else

		switch (corstr{idx}(1:4))
		case {'Beam'}
			[vel.inst(:,fdx,:)  ]  = VADCP.beam_to_instrument_STATIC(vel.beam(:,fdx,:),[]);
			[vel.ship(:,fdx,:)  ]  = ADCP.instrument_to_ship_STATIC(vel.inst(:,fdx,:),[],pitch(fdx),roll(fdx)-pi);
			[vel.earth(:,fdx,:) ]  = ADCP.ship_to_earth_STATIC(vel.ship(:,fdx,:),[],heading(fdx));
			disp('Beam coordinates found');
		case {'Inst'}
			[vel.ship(:,fdx,:)  ]  = ADCP.instrument_to_ship_STATIC(vel.inst(:,fdx,:),[],pitch(fdx),roll(fdx)-pi);
			[vel.earth(:,fdx,:) ]  = ADCP.ship_to_earth_STATIC(vel.ship(:,fdx,:),[],heading(fdx));
			disp('Instrument coordinates found');
		case {'Ship'}
			[vel.earth(:,fdx,:) ] = ADCP.ship_to_earth_STATIC(vel.ship(:,fdx,:),[],heading(fdx));
			disp('Ship coordinates found');
		case {'Eart'}
			disp('Earth coordinates found');
			% nothing to do
		otherwise
			error(['Unknown coordinate system ' corstr{idx}]);
		end % switch

		end % isempty btvel
	end % for idx
end % to_earth()

