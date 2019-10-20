% Tue Jul 29 08:43:28 WIB 2014
% Karl Kastner, Berlin
%% wrapper for conversion to beam velocity
%% Note that back-conversion to beam velocity is not unique in case of 3 beam
%% solutions (as RDI instruments doe not store which beams were used) and
%% if instrument internal bin-mapping is used (whichs precise algorithm remains
%% an RDI secret)
% TODO make beamangle optional argument
function [vel btvel] = to_beam(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad)
	if (isempty(btvel))
		btvel = zeros(length(FileNumber),4);
	end
	inverse = true;
	for idx=1:length(corstr)
		fdx = find(FileNumber == idx);
		switch (corstr{idx}(1:4))
			case {'Beam'}
				% nothing to do
			case {'Inst'}
				[vel.beam(:,fdx,:) btvel(fdx,:)] = HADCP.beam_to_instrument_STATIC(vel.inst(:,fdx,:),btvel(fdx,:),HADCP.DEFAULT_BEAMANGLE_RAD,inverse);
			case {'Ship'}
				[vel.inst(:,fdx,:) btvel(fdx,:)] = ADCP.instrument_to_ship_STATIC(vel.ship(:,fdx,:),btvel(fdx,:),pitch_rad(fdx),roll_rad(fdx),inverse);
				[vel.beam(:,fdx,:) btvel(fdx,:)] = HADCP.beam_to_instrument_STATIC(vel.inst(:,fdx,:),btvel(fdx,:),HADCP.DEFAULT_BEAMANGLE_RAD,inverse);
			case {'Eart'}
				[vel.ship(:,fdx,:) btvel(fdx,:)] = ADCP.ship_to_earth_STATIC(vel.earth(:,fdx,:),btvel(fdx,:),heading_rad(fdx),inverse);
				[vel.inst(:,fdx,:) btvel(fdx,:)] = ADCP.instrument_to_ship_STATIC(vel.ship(:,fdx,:),btvel(fdx,:),pitch_rad(fdx),roll_rad(fdx),inverse);
				[vel.beam(:,fdx,:) btvel(fdx,:)] = HADCP.beam_to_instrument_STATIC(vel.inst(:,fdx,:),btvel(fdx,:),HADCP.DEFAULT_BEAMANGLE_RAD,123,inverse);
			otherwise
				error(['Unknown coordinate system ' corstr{idx}]);
		end % switch
	end % for idx
end % to_beam()

