% Sat 30 Jun 19:06:53 CEST 2018
%% normalized particle radius
function [chi,k] = normalized_particle_radius(d_mm,f,varargin)
	d_m    = 1e-3*d_mm;

	% sound velocity water
	%c = 1500; % m/s
	if (~issym(d_mm))
		cw = Constant.sound_velocity_water(varargin{:});
	else
		syms cw
	end
	% wave length [lambda] = m
	% f = c/lamda
	lambda = cw./f;
	% wave number
	% [k] = 1/m
	k   = 2*pi./lambda;
	% normalized particle radius
	chi = k*d_m/2;
end
