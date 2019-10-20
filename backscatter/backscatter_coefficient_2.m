% Sat 30 Jun 19:06:53 CEST 2018
%
%% analytic basckatter coefficient
%% thorne 2002
%% thorne 2012
%
% function [ks2, fbs,chi] = backscatter_coefficient_2(d_mm,f,mode)
%
% [f_bs] = 1
% [ks^2] = [f^2/(r rho)] = m^2/kg	(thorne 2008, 1)
% [as]   = 1/m
% [d_mm] = mm
% f      = Hz = 1/s
%
% ks is also known as backscatter sensitivity
function [ks2, fbs, chi, D1] = backscatter_coefficient_2(d_mm,f,varargin)
	D    = 1e-3*d_mm;
	if (~issym(d_mm))
		rhos = Constant.density.quartz;
	else
		syms rhos
	end

	% critical diameter of transition to geometric range, chi = 1
	% D1 = 2/k;

	[fbs,chi] = backscatter_form_function(d_mm,f,varargin{:});
	% thorne 2008, 1
	% note that thorne-2014 has 3/(16 pi) in the system constant
	ks2 = 3/(16*pi)*fbs.^2./((D/2)*rhos);
end

