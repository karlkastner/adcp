% Tue 13 Jun 14:31:21 CEST 2017
% Karl Kastner, Berlin
%% wrapper for conversion to ship velocity
% TODO make beamangle optional argument
% TODO give btvel different fields
function [vel btvel] = to_ship_STATIC(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad)
	if (isempty(btvel))
		btvel = zeros(length(FileNumber),4);
	end
	inverse = true;
	for idx=1:length(corstr)
		fdx = find(FileNumber == idx);
		switch (corstr{idx}(1:4))
			case {'Beam'}
				[vel.inst(:,fdx,:) btvel(fdx,:)] = HADCP.beam_to_instrument_STATIC(vel.beam(:,fdx,:),btvel(fdx,:),HADCP.DEFAULT_BEAMANGLE_RAD,123,0);
				[vel.inst(:,fdx,:) btvel(fdx,:)] = HADCP.instrument_to_ship_STATIC(vel.ship(:,fdx,:),btvel(fdx,:),pitch_rad(fdx),roll_rad(fdx),0);
			case {'Inst'}
				[vel.inst(:,fdx,:) btvel(fdx,:)] = HADCP.instrument_to_ship_STATIC(vel.ship(:,fdx,:),btvel(fdx,:),pitch_rad(fdx),roll_rad(fdx),0);
			case {'Ship'}
				% nothing to do
			case {'Eart'}
				[vel.ship(:,fdx,:) btvel(fdx,:)] = HADCP.ship_to_earth_STATIC(vel.earth(:,fdx,:),btvel(fdx,:),heading_rad(fdx),inverse);
				[vel.inst(:,fdx,:) btvel(fdx,:)] = HADCP.instrument_to_ship_STATIC(vel.ship(:,fdx,:),btvel(fdx,:),pitch_rad(fdx),roll_rad(fdx),inverse);
			otherwise
				error(['Unknown coordinate system ' corstr{idx}]);
		end % switch
	end % for idx
end % to_ship_static

