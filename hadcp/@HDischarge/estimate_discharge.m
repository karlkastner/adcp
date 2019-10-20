% 2015-01-09 10:34:45 +0100
% Karl Kastner, Berlin
%% integrate and scale specifc discharge to total discharge for each ensemble
function [hq, hqerr, obj] = estimate_discharge(obj)

	obj.calc_specific_discharge_weights();

	% velocity to specific discharge
	c_u = obj.c_u;
	% specific discharge to total discharge
	c_q = obj.c_q;
	% scale up meausred velocity to total discharge
	w = c_u.*c_q;
	q_tot = w.*obj.velocity;
	% number of valid cells
	n = sum(isfinite(q_tot));
	% average
	hq    = nanmean(q_tot)';
	% standard error
	hqerr = (nanstd(q_tot)./sqrt(n-1))';
	
	if (0)
		% weighted average and standard error
		w = harmmean(w(~fdx))/n;
		hq = w*nansum(q_tot);
		hqerr = sqrt(w*nanvar(q_tot)/n);
	end

	% write back
	%obj.c_q   = c_q;
	obj.q_tot = q_tot;
	obj.Q     = hq;
	obj.Q_err = hqerr;
end % calc_total_discharge

