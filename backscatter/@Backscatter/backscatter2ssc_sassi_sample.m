% Sat 21 Jul 11:51:23 CEST 2018
%
%% convert backscatter to suspended sediment concentration
%% c.f. sassi
function  ssc = backscatter2ssc_sassi_sample(obj,bs0,Ibs0,bs_ref,Ibs_ref,p)
	if (nargin() < 6)
		p = obj.param;
	end
	ssc = (	bs0 ) ./ ( p(1).*bs_ref.^p(2) + p(3)*(Ibs0 - Ibs_ref) );
end
