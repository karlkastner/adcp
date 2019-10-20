% Thu Jun 27 16:45:53 UTC 2013
% Karl Kästner, Berlin
%%
%% acoustic cross sectin ? of sediment particles
%% Medwin, ch. 7.5.3
%% Axially Symmetric Spherical Mode Solutions
%
function [sigma_diff, sigma_total, L] = sigma_general(k,a_,theta,rho0,rho1,c0,c1)
	mmax = 1e4;
	mu  = cos(theta);
	sigma_diff  = zeros(size(a_));
	sigma_total = zeros(size(a_));
	L           = zeros(size(a_));
	for idx=1:length(a_)
	a = a_(idx);

	% number of terms to include
	m = ceil(k.*a+3);
	if (m > mmax)
		sigma_diff(idx)  = NaN;
		sigma_total(idx) = NaN;
		L(idx)           = NaN;
	else
		s_d = 0;
		s_t = 0;
		for mdx=0:m-1
			P_mdx = legendre_man(mu,mdx);
			%P_mdx = legendre(m,mu); P_mdx = P_mdx(end);
			C_mdx = C(k,a,rho0,rho1,c0,c1,mdx);
			% scattering cross section, medwin eq. 7.5.49
			% thorne 1992, eq 2
			s_d = s_d + (-1)^mdx*P_mdx*(2*mdx+1)/(1 + 1i*C_mdx);
			% total scattering cross section medwin eq. 7.5.51
			s_t = s_t + (2*mdx+1)/(1 + C_mdx^2);
		end
		L(idx)           = 1i*a./(k.*a).*s_d;
		sigma_diff(idx)  = abs(L(idx)).^2;
		sigma_total(idx) = 4*pi*a.^2./(k.*a).^2.*s_t;
	end
	end % for idx
end % function bs

function C_m = C(k,a,rho0,rho1,c0,c1,m)
	g  = rho1/rho0;
	h  = c1/c0;
	k1 = k/h;
	C_m =  ( jt(k1*a,m)*neumann_sphere(k*a,m) - g*h*bessel_sphere(k1*a,m)*nt(k*a,m) ) ...
            / ( jt(k1*a,m)*bessel_sphere(k*a,m) - g*h*bessel_sphere(k1*a,m)*jt(k*a,m) );
	function jt = jt(kR,m)
		jt = m/(2*m+1)*bessel_sphere(kR,m-1) - (m+1)/(2*m+1)*bessel_sphere(kR,m+1);
	end % function n
	function nt = nt(kR,m)
		nt = m/(2*m+1)*neumann_sphere(kR,m-1) - (m+1)/(2*m + 1)*neumann_sphere(kR,m+1);
	end % function nt
end % function C

