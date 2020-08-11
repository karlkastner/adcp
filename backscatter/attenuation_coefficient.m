% Sun  8 Jul 10:41:40 CEST 2018
%%
%% accoustic attenuation coefficient of suspended particles
%%
%% hanes 2012
%%
%% [d_mm] = mm
%% [f]    = Hz = 1/s
%% [as]   = 1/m (neper)
%% for db : chi_db  = 8.7 chi_neper
%% [M]    = kg/m^3 = mg/l
%% 
%% for normalization : chis = as(M=2650)
%% 
%% function [as,asnu,ass,X,chi] = attenuation_coefficient(d_mm,f,M,mode)
%%
function [as,asnu,ass,X,chi] = attenuation_coefficient(d_mm,f,M,mode,varargin)
	rhos =  Constant.density.quartz;
	if (nargin()<3 || isempty(M))
		% unit mass concentration (1kg/m^3)
		M = 1;
	end
	if (nargin()<4)
		%mode = 'thorne-theoretical';
		%mode = 'thorne-empirical';
		mode = 'moate-2012';
	end

	%[D] = m
	D    =  1e-3*d_mm;
	% [nu] = m/s
	nu   =  Constant.viscosity.kinematic.water;	
	rhow =  Constant.density.water;
%	c    =  Constant.sound_velocity_water(varargin{:});

	[chi,k] = normalized_particle_radius(d_mm,f,varargin{:});

	% Urick 1948, richards 2003, holdaway 1999, heathershaw 1996
	% viscous attenuation
	gamma = sqrt(pi*f/nu);
	s     = (9./(2*gamma*D)).*(1+2./(gamma*D)); 
	T     = 0.5 + 4.5./(gamma*D);
	asnu  =   ( 0.5*k.*M/rhos.*(rhos/rhow-1)^2 ...
		.*( s./(s.^2+(rhos/rhow+T).^2) ) );
	
	% attenuation by scattering (shadowing)
	switch (mode)
	case {'urick-1948'}	
		% ass = (D*M*k^4)/(12*pi*rhos);
		X   = 1/4.5*chi.^4;
	case {'thorne-manley-1993'}
		if (0)
			Ka = (gk^2 + gr^2/3)/6;
		else
			% from thorne-hardcastle 1993
			Ka = 0.18;
		end
		X  = (4/3)*Ka*chi.^4./(1 + 0.9*chi.^(5/2) + 4/3*Ka*chi.^4);
	case {'thorne-hardcastle-1993'}
		ka = 0.18;
		X = 4/3*ka*chi.^4./(1 + chi.^2 + 4/3*ka*chi.^4);
	case {'medwin-clay-1997'} % c.f. RUSSO and BOSS 2012
		% scattering cross section
		sigma = TODO;
		X = 4/pi*sigma;
	case {'thorne-hanes-2002'}
		ka = 0.18;
		X = 1.1*(4/3)*ka*chi.^4./(1 + 1.3*chi.^2 + 4/3*ka*chi.^4);
	case {'thorne-meral-2008'}
		% thorne 2008 eq 9
		% hanes 2012
		% scattering cross-section
		X   = ( 0.29*chi.^4 )./(0.95 + 1.28*chi.^2 + 0.25*chi.^4);
		%ass = 1.5./D.*M/rhos.*X;
	case {'betterridge-thorne-2008'}
		X = 0.24*(1-0.4*exp(-((chi-5.5)/2.5).^2)) ...
			.*chi.^4./(0.7 + 0.3*chi + 2.1*chi.^2 - 0.7*chi.^3 + 0.3*chi.^4);
	case {'thorne-theoretical'}
		% thorne 2008, eq 5b
		e  = 39;
		%  X ~ 0.26
		X  = 2*( ((e-1)/(3*e))^2 + 1/3*((rhos - rhow)/(2*rhos+rhow))^2);
		X  = X*chi.^4;
		%f   = 1.17*chi^2;
		%ass = 1.5./D.*M/rhos.*X;
	case {'moate-2009'}
		% asymptote for raleigh   0.32 xi^4
		% asymptote for geometric 1.4
		X = 0.29*chi.^4./(0.92 + 0.9*chi.^2 + 0.21*chi.^4);
	case {'moate-2012'} % c.f. thorne 2014, aberle 2017
		X = rhos*(0.09*chi.^4)./(1380 + 560*chi.^2 + 150*chi.^4);
	otherwise
		error('here');
	end
	% thorne 2008, 1, a3
	% 1.5 neper ~ 6.5 db, c.f guerrero 2011
	ass = 3/2*M./(D*rhos).*X;

	as = asnu + ass;
end

