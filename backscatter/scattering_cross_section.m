% 2020-07-26 11:16:51.827004804 +0800 adcp/backscatter/scattering_cross_section.m
% TODO, this is just part of the backscatter form function


thorne mereal
	% a : radius
	x  = 2*pi*a./f;
	xi = 0.29.*x.^2 ./ ( 0.95 + 1.28*x.^2 + 0.25*x.^4);
