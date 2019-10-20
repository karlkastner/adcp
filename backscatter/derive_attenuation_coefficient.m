% 2018-07-08 17:05:32.226840068 +0800
syms D M rhos k xi

		v = 4/3*pi*(D/2)^3
		n = M/(rhos*v)
		% Urick (4)
		ass = 2/9*k^4*(D/2)^4*pi*(D/2)^2*n
xi_ = k*D/2
		ass = subs(ass,k,solve(xi-xi_,k))
