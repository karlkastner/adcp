% Thu Jun 27 16:55:37 UTC 2013
% Karl Kästner, Berlin

%% differential cross section
%% geometrical backscattering for spherical bodies
%% ka >> 1, large particles or high frequencies
%
%% k : wave number
%% a : radius of the particle
function [sigma_diff, sigma_total] = bs_geometric(k,a)
	if (k*a < 1)
		sigma_diff = NaN;
		sigma_total = NaN;
	else
		sigma_diff = 0.25*a^2;
		% the total geometric scattering cross section
		% equals the area of the cross section
		sigma_total= pi*a^2;
	end
end % bs_geometrical scattering

