% 2014-10-23 11:07:05.248813053 +0200
%% rotate ensembles in stream direction (transverse velocity integrates to zero)
function [u_, v_] = streamvel(u,v)
	mu		 = nanmean(u);
	mv		 = nanmean(v);
	hyp		 = hypot(mu,mv); %sqrt(mu*mu + mv*mv);
	cos_		 = mu./hyp;
	sin_		 = mv./hyp;
	u_		 = cos_*u + sin_*v;
	v_		 =-sin_*u + cos_*v;	
end

