% Mon 17 Jul 16:05:14 CEST 2017
%
%% convert acoustic backscatter to sediment concentration
function [C C_ beta_] = backscatter_to_concentration(Z,beta)
	a   = 0.17;
	b   = 0.24; 
%	chi = 1.2;
	chi = 120;
%	chi = 400;
%	chi = 30;
%	chi = 60;
%	chi = 90;
%	chi = 105;
%nanmedian(beta(:))
%pause
	r = 1.0;
	beta_r = interp1(Z,beta,r,'linear');
	K = a*beta_r.^b;
	gamma = log(10)/5*chi;
	% TODO, true integral (blanking range)
	% TODO integrate from r, not 0
	dZ = Z(2)-Z(1);
	Ibeta = cumsum(beta)*dZ;
	% the denominator is ill conditioned and should be regularised
	den = bsxfun(@minus,K,gamma*Ibeta);
	den = medfilt1(double(den)',10)';
	
	C = beta./max(0,den);
	% bsxfun(@minus,K,gamma*Ibeta));
	% simple calibration, no attenuation
	C_ = bsxfun(@times,beta,1./K);

	beta_ = bsxfun(@times,C,1./K);
end

