% Thu Jun 27 16:48:03 UTC 2013
% Karl Kästner, Berlin
%
%% Rayleigh scattering for a sphere (ka << 1)
%% small particles or low frequencies
%% Medwin 7.5.2 Rayleigh Scatter From a Sphere (ka << 1)
function [sigma_diff sigma_total] = bs_rayleigh(k,a,theta)
	if (k*a > 1)
		sigma_diff  = NaN;
		sigma_total = NaN;
	else
		sigma_diff  = 1/9*(k*a)^4*a^2*(1 + 1.5*cos(theta))^2;
		sigma_total = a^2 * 25/36 * (k*a)^4;
	end
end % function bs_rayleigh

