% 2018-06-07 19:20:10.393760237 +0200
%% convert backscatter to suspended sediment concentration, implicit method
%
function [c,lbs_] = backscatter2ssc_sample_implicit(obj,R,S,lbs,S0,last,param)
	if (nargin()<7)
		param = obj.param;
	end
	for idx=1:size(lbs,2)
		% determine attenuation corrected backscatter and
		% concentration
		[ci,lbsi]   = backscatter2ssc_implicit(R(1:last(idx),idx),lbs(1:last(idx),idx),param);
		% interpolate to sample depth
		c(idx,1)    = interp1(S(1:last(idx),idx),ci,S0(idx),'linear','extrap');
		lbs_(idx,1) = interp1(S(1:last(idx),idx),lbsi,S0(idx),'linear','extrap');
	end % for idx
end % bs2ssc

