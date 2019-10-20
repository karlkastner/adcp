% Sun Mar  2 16:07:55 WIB 2014
% Karl Kastner, Berlin
%% z-mapping, i.e. correct for roll and pitch of instrument
function [bs, dir_earth] = map_z(bs, R, beamangle, pitch, roll)

	% directions of the beams in instrument coordinates
	% columns are beam directions [dx,dy,dz]'
	s = sin(beamangle);
	c = cos(beamangle);
	% TODO, function for this matrix
	dir_inst = [s -s 0  0
	            0  0 s -s
	            c  c c  c];

	% directions of the beams in earth coordinates (heading ignored)
	% inverse rotation matrix for roll
	%roll = median(roll);
	cr = cos(roll);
	sr = sin(roll);
	RRi = [cr 0 sr;
		0 1  0;
	      -sr 0 cr];

	% inverse of rotation matrix for pitch
	%pitch = median(pitch);
	pitch = atan(tan(pitch).*cos(roll));
	cp    = cos(pitch);
	sp    = sin(pitch);
	RPi   = [1  0   0;
	         0 cp -sp;
		 0 sp  cp];
	RR = RPi * RRi;
	%RR = inv(RR);
	dir_earth = RR*dir_inst;
	
	% backscatter vs true z-value of bins
	% interpolate backscatter to equally distributed z-values
	if (2 == length(size(bs)))
		for idx=1:size(bs,2)
			bs(:,idx) = interp1(R*dir_earth(3,idx),bs(:,idx),R,'linear');
		end
	else
		for idx=1:4
			bs(:,:,idx) = interp1(R*dir_earth(3,idx),bs(:,:,idx),R,'linear');
		end
	end

end % map_z

