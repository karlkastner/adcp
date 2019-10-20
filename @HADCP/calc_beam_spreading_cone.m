% Sun Jun  8 14:55:51 WIB 2014
% Karl Kastner, Berlin
%
%% beam spreading
%% Note: beams spread in the form of bessel functions
%%       this is the engineering approach as cones, which is however not
%%       a good approximation, it is better to approximate it as a gaussian
function obj = calc_beam_spreading_cone(obj)
	theta = obj.beamangle_rad;
	d = [ cos(theta) 1  cos(theta)
	      sin(theta) 0 -sin(theta)
	      0          0           0];
	% scale side beams
	s = sqrt(1+sin(theta)^2);
	d(:,1) = s*d(:,1);
	d(:,3) = s*d(:,3);
	
	for idx=1:length(pitch)
		p = obj.pitch(idx);
		r = obj.roll(idx);
		% sensor pitch to true pitch
		p = atan(tan(p)*cos(r));
		P = [ cos(p) 0  -sin(p)
		      0      1       0
		      sin(p) 0  cos(p)];
		R = [ 1      0       0;
		      0  cos(r) sin(r)
		      0 -sin(r) cos(r) ];
		d_ = P*R*d;
	
		% beam centre difference
		Dz(idx,:) = d_(3,:);
		% beam top difference
		% beam bottom difference
	end % for idx
end
	
