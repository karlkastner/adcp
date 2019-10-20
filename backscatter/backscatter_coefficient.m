% % Thu Jul  4 11:22:59 UTC 2013
%% analytic determination of the backscatter coefficient
function [ks2,xi] = backscatter_coefficient(d,f)

	% sound velocity in water
	c0   = Constant.sound_velocity_water();

	% sound velocity in quartz
	c1   = Constant.sound_velocity_quartz;

	% density of water
	rho0 = Constant.density.water;

	% density of quartz
	rho1 = Constant.density.quartz;

	% scattering angle
	theta = 0;

	% wave length in water
	lambda = c0/f; % m

	% wave number in water
	k=2*pi/lambda;

	% radius
	a = 0.5*d;

	% Sassi, eq 2
	% mass concentration of particles
	% M = V*n
	% S = 10 log_10(ks2*M)
	% S = 10 log_10(n*sigma)
	% => ks2 = sigma/V
	%    M = 1/ks2*10^(S/10)
	
	% scattering cross section
	[ds, ts] = scattering_cross_section_general(k,a,theta,rho0,rho1,c0,c1);

	% target strength ?
	tsn = ts./(pi*a.^2);

	% volume of sphere
	V = 4/3*pi*a.^3;
	
	% attenuation coefficient ?
	% thorne 1992, eq 11
	xi = (3./(4*a*rho1)).*((a.*(a.^2.*tsn))./(a.^3));

	% squared form function
	% eq 3 in Sassi 2013
	% ds = 1/4*f^2
	% this was ds
	f2 = 4*ts./(a.^2);

	% squared sediment constant (suspension parameter)
	% eq 5 in Sassi : TODO these equations differ
	%ks2 = ds./V;
	%ks2 = (a.*(a.^2.*f2)./(a.^3))./(a*rho1)
	%ks2 = 3./(16*pi*rho1).*(a.^2.*f.^2)./(a.^3)
	ks2 =  

end

