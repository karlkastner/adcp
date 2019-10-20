% Mon Jan 27 13:52:41 WIB 2014
% Karl Kastner, Berlin

%% correct the bottom coordinates for pitch and roll
function obj = correct_coordinates(obj)
	% ship coordinates (origin of beams)
	C0 = [X, Y, 0];

	beamangle = obj.beamangle_rad;

	% length of the beam before hitting the bottom
	L = (1/cos(deg2rad(beamangle_rad)))*Z;
	for idx=1:4
		% H*(-P)*(-R)*[1 0 0]

		% direction of the beam
		dir = [   ];
		% location of beam hitting the ground
		C(:,:,idx) = repmat(L(:,idx),1,3).*dir + C0;
	end

	% write values
	%obj.C = C;
end % correct_coordinates

