% Sat 30 Jun 19:06:53 CEST 2018
%
%% acoustic backscatter form function
%
function [fbs,chi] = backscatter_form_function(d_mm,f,mode,varargin)
	if (nargin()<3 || isempty(mode))
		%mode = 'thorne-2014';
		mode = 'thorne-meral-2008';
		%mode = 'medwin-clay-1977';
	end

	chi = normalized_particle_radius(d_mm,f,varargin{:});

	if (~issym(d_mm))
	rhos =  Constant.density.quartz;
	rhow =  Constant.density.water;
	else
		syms rhos rhow
	end

	switch (mode)
	case {'johnson-1977'}
		% fbs = sigms/(pi a^2) ?
		%alpha = 4/9*chi^4*(1 - b + 3*(g-1)/(1+2*g))^2
		%g = 
		%h = 
		alpha  = 4/9*((1-g*h^2)/(3*g*h^2) + (1-g)/(1+2*g))^2
		fbs2   = 2*alpha*chi.^4./(2 + 3*chi.^4);
		fbs    = sqrt(fbs2);
	case {'medwin-clay-1997'}
		% c.f. thorne 2008
		e   = 39;
		% f ~ 1.17
		fbs_ = 2*( (e-1)./(3*e) + (rhos-rhow)./(2*rhos+rhow) );
		fbs  = fbs_*chi.^2;
	case {'hay-sheng-1988'}
		% kf = 2/3*(gk - gr); % for cos theta = 1
		kf = 1.1; % from thorne 1993
		fbs = kf*chi.^2./(1 + kf*chi.^2);
 	case {'thorne-manley-1993'} % c.f. russo and boss 2012
		% theoretic equation
		nu1  = 0.37;
		eta1 = 0.5;
		x1 = 1.4;
		nu2 = 0.28;
		eta2 = 2.2;
		x2 = 2.8;
		c0 =   (1-nu1*exp(-((chi-x1)/eta1).^2)) ...
		     .*(1+nu2*exp(-((chi-x2)/eta2).^2));
		% compressibility (inverse of bulk modulus), values not given in reference
		if (0)
		kw  = 1/2.2;
		ks  = 1/38;
		gk  = (ks-kw)/kw;
		gr  = 3*(rhos-rhow)/(2*rhos+rhow);
		kf  = 2/3*(gk - gr);
		else
			kf = 1.1; % from thorne 1993
		end
		fbs = c0.*(kf*chi.^2)./(1+kf*chi.^2);
	case {'thorne-hardcastle-1993'}
		nu1  = 0.25;
		x1 = 1.4;
		eta1 = 0.5;
		nu2 = 0.37;
		eta2 = 2.2;
		x2   = 2.8;
		kf = 1.1;
		c0 =   (1-nu1*exp(-((chi-x1)/eta1).^2)) ...
		     .*(1+nu2*exp(-((chi-x2)/eta2).^2));
		fbs = c0.*(kf*chi.^2)./(1+kf*chi.^2);
	case {'thorne-hanes-2002'}
	fbs = 1.1*(1 - 0.25*exp( -(2*(chi-1.4)).^2 )) ...
		.*(1 + 0.37*exp( -((chi-2.8)/2.2).^2 )) ...
		.*(1.1*chi.^2./(1+1.1*chi.^2));
	case {'thorne-meral-2008'}
	% nb : tm give also a theoretic equation of the scattering
	% raleig asymptote: f=1.17x^2
	% geometric asymptote: f=1.1 (const)
	fbs =       (1-0.35*exp(-((chi-1.5)/0.7).^2)) ...
		  .*(1+0.50*exp(-((chi-1.8)/2.2).^2)) ...
		  .*chi.^2./(1 + 0.9*chi.^2);
	case {'moate-2009'}
		% raleigh asymptote: fbs = 1.26 xi^2
		% geometric asymptote: fbs = 1.8 (const)
		fbs = (1-0.21*exp(-((chi-1.6)/0.55).^2)) ...
		      .*(1+0.30*exp(-((chi-0.3)/0.8).^2)) ...
		      .*(1+0.53*exp(-((chi-2.55)/0.95).^2)) ...
			.*chi.^2./(1+0.56*chi.^2);
	case {'hurther-2011'}
		% there seems to be a mistake in the reference
		D = 1e3*d_mm;
		fbs = 1.1*1*(1.1*chi.^2/(1 + 1.1*chi.^2))*D;
	case {'betterige-thorne-2008'}
		fbs =    (1-0.5*exp(-((chi-1.5)/0.5).^2)) ...
		       .*(1+0.4*exp(-((chi-1.5)/3.0).^2)) ...
		       .*(1-0.5*exp(-((chi-5.9)/0.7).^2)) ...
			.*chi.^2./(1.17 + 0.95*chi.^2);
	case {'moate-thorne-2012'}
		% they give two variants, why?
		if (0)
		%c1 = 0.24; c2=1.19; c3=0.56; c4=0.13;
		d1 = 0.75;
		d2 = 1.4;
		d3 = 1.35;
		d4 = 2.1;
		d5 = 1.45;
		d6 = 1.25;
		d7 = 0;
		d8 = 0;
		d9 = 1;
		d10 = 1;
		d11 = 0.41;
		fbs =      (1-d1*exp(-((chi-d2)/d3).^2)) ...
			.* (1+d4*exp(-((chi-d5)/d6).^2) ...
			.* (1+d7*exp(chi-d8)/d9).^2) ...
			.*chi.^2 ./ (d10 + d11*chi.^2);
		else
			fbs = sqrt(rhos).*(1-0.25*exp(-((chi-1.5)/0.35).^2)) ...
			         .*(1+0.60*exp(-((chi-2.9)/1.15).^2)) ...
				 .*chi.^2./(42 + 28*chi.^2);
		end
	case {'thorne-2014'} % TODO, there is a mistake in the source, 
			     % apparently thorne changed the notation somewhat K != ks
	fbs = sqrt(rhos) ...
               * (1 - 0.25*exp(-((chi-1.5)/0.35).^2 )) ...
              .* (1 + 0.6*exp(-((chi-2.9)/1.15).^2)) ...
	      .* chi.^2./(42 + 25*chi.^2);
	otherwise
		error('here');
	end
end % backscatter_form_function

