% 2015-03-11 10:08:52.449062265 +0100
%% homogenise the horizontal velocity profile
function [scale, U_scaled] = homogenise_profile(val,m2)

	% average over time (all velocity profiles)
	if (nargin()<2)
		m2 = nanmean(val,2);
	end

	% scales, inverse of the joint velocity profile
	scale = nanmean(m2)./m2;

	% scaled velocity
	if (nargou()>1)
		Us = bsxfun(@times,scale,U_scaled);
	end
end % homogenise_profile

