% Mon Sep 30 10:16:31 UTC 2013
% Karl Kästner, Berlin
%
%% backscatter inversion
%%
function obj = bsinvert(obj, jdx)
	nt = size(obj.beta_f,2);
	nr = length(obj.R);
	%beta_k_f = repmat(obj.beta_f(obj.nk,:),nr,1);
	%I_k_f = repmat(obj.I_f(obj.nk,:),nr,1);
	beta_k_f = repmat( interp1(obj.D,obj.beta_f,obj.d_k),nr,1);
	I_k_f = repmat( interp1(obj.D,obj.I_f,obj.d_k),nr,1);
	% convert backscatter to sediment concentration
	if (0 == jdx)
		pmean = obj.param.mean;
		pcov  = obj.param.cov;
	else
		pmean = obj.J(jdx).param.mean;
		pcov  = obj.J(jdx).param.cov;
	end

	obj.M.mean = feval(obj.funcC{obj.model}, pmean, I_k_f, obj.I_f, beta_k_f, obj.beta_f, [], []); %obj.PC, obj.D);

	sb = repmat(obj.sbeta,1,nt);
	sB = repmat(obj.sI,1,nt);
	B = obj.I_f - I_k_f;
	M = obj.M.mean;
	switch (obj.model)
		case {2}
			a = pmean(1);
			b = pmean(2);
			beta = obj.beta_f;
			obj.M.bias = 0.5*M.*( ...
				  log(beta).^2*pcov(1,1) ...
				+ b*(b-1)./beta.^2.*sb.^2 ...
				+ 2*log(beta)/a*pcov(1,2) );
			obj.M.std  = M.*sqrt( ...
                                  1/a^2*pcov(1,1)  ...
				+ b^2./beta.^2.*sb.^2 ... 
				+ log(beta).^2*pcov(2,2) ...
                                + 2*log(beta)/a*pcov(1,2) );
		case {5}
			% at the moment, there is no good estimate for sB, so assumed to be zero
			obj.M.bias = obj.M.mean.^3./obj.beta_f.^2.*( ...
					pcov(1,1) ...
					- 2*pcov(1,2)*B ...
					+ pcov(2,2)*B.^2 ...
					+ sB.^2*pmean(2)^2);
			obj.M.std  = sqrt( (obj.M.mean./obj.beta_f).^2.*sb.^2 ...
                	           + (obj.M.mean.^4./obj.beta_f.^2).*(pcov(1,1) - 2*pcov(1,2)*B + pcov(2,2)*B.^2 + sB.^2*pmean(2)^2) );
		case {11}
			B = obj.I_f;
			obj.M.bias = 0.5*M.*(pcov(1,1) + pcov(2,2)*B.^2 + 2*pcov(1,2)*B);
			obj.M.std  = M.*sqrt(pcov(1,1) + pcov(2,2)*B.^2 + 2*pcov(1,2)*B);
		otherwise
			obj.M.bias = zeros(size(M));
			obj.M.std  = zeros(size(M));
	end
end % invert()

%
% auxilary functions
%

% determine standard error of the reconstructed sediment concentration
% with respect to changes of the calibration points
function [M_c, sM] = sfunc(func, param, Jparam1, Jparam2, varargin)
	n = size(Jparam1,1);
	n2 = size(Jparam2,1);

	% get the estimate with parameters based on all calibration points
	M   = func(param, varargin{:});

	% get the all estimates leaving out one calibration point
	M_J1 = [];
	for idx=1:n
		M_J1(:,:,idx) = func(Jparam1(idx,:), varargin{:});
	end

	if (isempty(Jparam2))
		% get a first order bias corrected estimate
		% theta_c =  n*theta - (n-1)/n sum theta_J
		%M = log(M); M_J1 = log(M_J1);
		M_c = n*M - (n-1)/n*sum(M_J1,3);
		%M_c = exp(M_c);
	else
		% get all the estimates leaving out two calibration points
		M_J2 = [];
		for idx=1:n2
			M_J2(:,:,idx) = func(Jparam2(idx,:), varargin{:});
		end
		% get a second order bias corrected estimate
		%M = log(M); M_J1 = log(M_J1); M_J2 = log(M_J2);
		M_c = 1/(2*n-1)*(n^3*M - (2*n^2-2*n+1)*(n-1)/n*sum(M_J1,3) + (n-1)^2*(n-2)*2/(n*(n-1))*sum(M_J2,3));
		%M_c = exp(M_c);
		%M_c = 1/(n-1)*(n^3*M - (2*n^2-2*n+1)*(n-1)/n*sum(M_J1,3) + (n-1)^2*(n-2)*2/(n*(n-1))*sum(M_J2,3));
	end

	% estimate the variance
	% TODO, this is not the bias corrected variance
	sM = std(M_J1,[],3);
end % sfunc

