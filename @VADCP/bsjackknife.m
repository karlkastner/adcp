% Wed Oct 16 11:36:07 UTC 2013
% % Karl Kästner, Berlin
%% compute the jackknife estimates of the parameters and their covariances
%
function bsjackknife(obj,J)
	nc = length(obj.M_ref);
	np = length(obj.param.mean);
	for idx=1:length(J)
		nj = J(idx);
		% get all subset indices with nj samples deleted
		P = nchoosek(1:nc,nc-nj);
		m = size(P,1);
		obj.J(nj).P = P;
		obj.J(nj).param.mat = zeros(size(P,1),np);
		% compute all subset parameter sets with nj samples deleted
		for jdx=1:m
			pdx = P(jdx,:);
			[pmean M_reg err A] = feval(obj.regC{obj.model}, ...
				obj.I_f_1(pdx), ...
				obj.I_f_2(pdx), ...
				obj.beta_f_1(pdx), ...
				obj.beta_f_2(pdx), ...
				obj.M_ref(pdx), ...
				ones(nc,1)*obj.d_k, ...
				obj.r_ref(pdx), ...
				obj.PC_ref);
			obj.J(nj).param.mat(jdx,:) = pmean(:)';
		end
		% determine the pseudo-parameter
		% TODO, this is only valid for leave-out one
		n = nc;
		d = nj;
		obj.J(nj).param.pseudo = n/d*repmat(obj.param.mean,m,1) - (n-d)/d*obj.J(nj).param.mat;
			%obj.J(nj).param.pseudo = m*repmat(obj.param.mean,m,1) - (m-1)*obj.J(nj).param.mat;
		% get the jackknife estimate of the mean
		obj.J(nj).param.mean = n/d*obj.param.mean - (n-d)/d*mean(obj.J(nj).param.mat);
			%obj.J(nj).param.mean = mean(obj.J(nj).param.pseudo);
			%obj.J(nj).param.mean = m*obj.param.mean + (m-1)*mean(obj.J(nj).param.mat);
		% get the jackknife estimate of the covariance matrix
			% obj.J(nj).param.cov  = 1/m*cov(obj.J(nj).param.pseudo);
		%obj.J(nj).param.cov = (n-d)/(d*nchoosek(n,d))*cov(obj.J(nj).param.mat);
		obj.J(nj).param.cov = (n-d)/(d)*cov(obj.J(nj).param.mat);
	end
end % bsjackknife

