% Sat 21 Jul 14:40:39 CEST 2018
%
%% convert suspended sediment concentration to backscatter
% 
%% function bs = ssc2backscatter(ssc,d_mm,f,varargin)
%%
%% ssc : mass concentration of sediment [ssc] = g/l = kg/m^3
%% d_mm : grain size diameter [d_mm] = mm
%% f : frequency [f] = Hz = 1/2
function bs = ssc2backscatter(ssc,d_mm,freq,varargin)
	rhos = Constant.density.quartz;
	D = 1e-3*d_mm;
	% mass of a single particle
	m1 = rhos*Geometry.volume_sphere(D/2);
	% number of particles
	n = ssc/m1;
	% form function
	f = backscatter_form_function(d_mm,freq,varargin{:});



	% scattering cross section
	sigma = 1/4*(D/2)^2*f^2;
	% backscatter
	% bs = ks^2 ssc
	bs = n*sigma;

	% note that this makes it consistend with Thorne meral, eq 1
	% sassi 1212 has the spurious factor of 3/(16 pi) in 
	bs = bs/(3/(pi*16));
end

