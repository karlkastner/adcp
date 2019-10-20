% Mon Sep 30 10:17:09 UTC 2013
% Karl Kästner, Berlin
%
%% backscatter to sediment calibration
%%
%% calibtation subroutine
%% M_ref : sediment concentration calibration values 
%% d_k   : depth of virtual reference value K
%%         (choose close to receiver, but out of near field, e.g. within 2m .. 4m)
%% TODO : better documentation of input values
%% TODO : rename nk into ik, bacause it is an index and not a length
%% TODO rename r_ref and d_k into r_1 and r_2
function obj = bscalibrate(obj, model, M_ref, d_ref, t_ref, d_k)
	% get array lengths
	nc = length(M_ref);
%	nr = length(obj.R);
%	nt = length(obj.temp);

	% take over input variables
	obj.model   = model;
	obj.M_ref   = M_ref;
	obj.d_ref   = d_ref;
	obj.t_ref   = t_ref;
	obj.d_k     = d_k;
	obj.model   = model;

	if (0)
	% convert depth to radial distance of calibration points
	obj.r_ref = d_ref ./ cos(obj.beamangle);
	
	% find space and time indices of reference values
	obj.nr_ref = zeros(length(d_ref),1);
	obj.nt_ref = zeros(length(d_ref),1);
	for idx=1:length(d_ref)
		%obj.nr_ref(idx,1) = find( (obj.D - d_ref(idx)) > 0,1);
		obj.nr_ref(idx,1) = round( max((d_ref(idx) - obj.L),0) / obj.dr );
		obj.nt_ref(idx,1) = find( (obj.time - t_ref(idx)) > 0,1);
	end
	% depth of virtual reference value, this value is fixed
	obj.nk = max(1, round( max((d_k - obj.L),0) / obj.dr ));

	% flat indices of regression parameters
	i1dx = sub2ind([nr nt], obj.nk*ones(size(obj.nt_ref)), obj.nt_ref);
	i2dx = sub2ind([nr nt], obj.nr_ref, obj.nt_ref);
	% 1 : virtual reference point
	% 2 : actual calibration point
	obj.I_f_1    = obj.I_f(i1dx);
	obj.I_f_2    = obj.I_f(i2dx);
	obj.beta_f_1 = obj.beta_f(i1dx);
	obj.beta_f_2 = obj.beta_f(i2dx);
	% obj.beta_k_f = repmat(obj.beta_f(obj.nk,:),nr,1);
	% obj.I_k_f = repmat(obj.I_f(obj.nk,:),nr,1);

	else
		obj.I_f_1 = interp2(obj.D, obj.time, obj.I_f', repmat(d_k, nc,1), t_ref,'linear');
		obj.I_f_2 = interp2(obj.D, obj.time, obj.I_f', obj.d_ref, t_ref,'linear');
		obj.beta_f_1 = interp2(obj.D, obj.time, obj.beta_f', repmat(d_k, nc,1), t_ref,'linear');
		obj.beta_f_2 = interp2(obj.D, obj.time, obj.beta_f', obj.d_ref, t_ref,'linear');
			% this is for plotting access internally interpolation should be used
			% obj.nr_ref = round( interp1(obj.D,(1:length(obj.D))',obj.d_ref, 'linear') );	
	end

%	if (13 == model)
%		I_k_f = repmat(obj.I_f(obj.nk,:),nr,1);
%		B = obj.I_f - I_k_f ;
%		m = 3;
%		tic()
%		[obj.PC obj.g_K A AB] = prepare13(m, obj.R, obj.beta_f,B);
%		% check convexity
%		x = obj.g_K*linspace(0.5,2,41);
%		x = logspace(-4,4,40);
%		for idx=1:length(x)
%			r = fopt13(x(idx), A, AB, obj.beta_f(9:end-8,:));
%			nr(idx) = norm(r);
%		end
%		figure(1010); clf();
%		%semilogx(x,nr)
%		plot(x,nr);
%		hold on;
%		plot(obj.g_K,1e-3,'*');
%		hold off;
%		obj.g_K
%		toc()
%		pause
%		[] = obj.PC(:,obj.nt_ref)';
%	else
%		[] = zeros(size(obj.r_ref));
%	end


	% compute calibration parameters with least square regression
	% using the respective backscatter model

	%[obj.param obj.M_reg obj.err] = regfunc([I_1 I_2 beta_1 beta_2 M_2 r_1*ones(n,1) r_2] );
	%[obj.M_ref obj.beta_f_1 obj.beta_f_2 obj.I_f_1 obj.I_f_2]
	try
	[obj.param.mean, obj.M_reg, obj.err A] = feval(obj.regC{model}, obj.I_f_1, obj.I_f_2, obj.beta_f_1, obj.beta_f_2, obj.M_ref, ones(nc,1)*d_k, [], []);
	catch e
		disp(e);
		obj.rmse = NaN;
		% TODO, set other values as well to 0
		return;
	end
	obj.param.mean = obj.param.mean(:)';

%	% extract concentration values at calibration points
%	obj.M_reg = obj.M(idx);
%	% error of the regression
%	obj.err = obj.M_reg - obj.M_ref;

	% Pearsson correlation (alternative Spearman)
	obj.corrMM = corr(obj.M_ref, obj.M_reg);
	obj.corred = corr(obj.err, obj.d_ref);
	obj.corret = corr(obj.err, obj.t_ref);

	% L2-norm of the error at the calibration points
	obj.nerr = sqrt(obj.err'*obj.err);
	
	% bias at the calibration points
	% for a linear regresion, this has to approach 0
%	obj.bias = mean(obj.err);
%	obj.bias = exp(mean(log(obj.M_ref) - log(obj.M_reg)));
%	obj.bias = mean(log(obj.M_ref) - log(obj.M_reg));

	% root mean square error at the calibration points 
	obj.rmse = sqrt(1/nc*obj.err'*obj.err);

	% chi^2-test for goodness of fit
	chi2  = obj.err'*obj.err/var(obj.M_ref);
	np  = length(obj.param);
	obj.chi2r = 1/(nc-np-1)*chi2;
	% normality of errors
	try
		obj.abnormal  = ttest(obj.err);
		%obj.abnormal = chi2gof(obj.err);
		%obj.abnormal = kstest( (obj.err - mean(obj.err))/std(obj.err) )
		% lillietest(obj.err)
	catch
		obj.abnormal = NaN;
	end

	% jacknife leave-one out matrix, each row contains the set of calibration
%	obj.Jparam  = jackknife(obj.reg, obj.I_f(i1dx), obj.I_f(i2dx), obj.beta_f(i1dx), obj.beta_f(i2dx), M_ref, ones(nc,1)*d_k, obj.r_ref, []);
%	obj.Jparam  = jacknife2(obj.reg, obj.I_f(i1dx), obj.I_f(i2dx), obj.beta_f(i1dx), obj.beta_f(i2dx), M_ref, ones(nc,1)*d_k, obj.r_ref, []);

	% jacknife leave-two out matrix
	%obj.Jparam2 = jackknife2(regfunc, I_1, I_2, beta_1, beta_2, M_2, ones(n,1)*r_1, r_2);
	%obj.Jparam2  = jackknife2(obj.reg, obj.I_f(i1dx), obj.I_f(i2dx), obj.beta_f(i1dx), obj.beta_f(i2dx), M_ref, ones(nc,1)*d_k, obj.r_ref);

	% get an estimate of the regression parameters
	% with a first order correction of the _statistical_ bias
        % theta_c =  n*theta - (n-1)/n sum theta_J                      
%	obj.jparam.mean = nc*obj.param.mean - (nc-1)/nc*sum(obj.Jparam,1);

	% variance-covariance matrix of the regression parameters
	% (diagonal is square of standard error)
	% without jackknife
	np = length(obj.param.mean);
	s2 = 1/(nc-np)*obj.err'*obj.err;
	obj.param.cov = s2*inv(A'*A);
	obj.param.std = sqrt(diag(obj.param.cov));

	% 2-sigma (95%) confidence intervals of the parameters
	% TODO nc-1 or nc-np-1 ?
%	s = diag(obj.sparam);
%	obj.cparam = [obj.param(:) - tinv(0.975,nc-1)*s(:), ...
%                     obj.param(:) + tinv(0.975,nc-1)*s(:)];

	% correlation of the parameters with respect to changes
	% in the calibration points
	% obj.param.corr = obj.param.cov(1,2)./sqrt(obj.param.cov(1,1)*obj.param.cov(2,2));
	D = diag(1./sqrt(diag(obj.param.cov)));
	obj.param.corr = D*obj.param.cov*D;
end % calibrate()

% r or d ?
function [x g_K A AB] = prepare13(m, r, beta, B)
	% exclude the first and last 4m
	r = r(9:end-8);
	beta = beta(9:end-8,:);
	B = B(9:end-8,:);

	nr = size(beta,1);
	nt = size(beta,2);
%	PC = zeros(2*m,nt);

	% set up the constant part of the regression matrix
	A = zeros(nr,m);
	A(:,1) = ones(nr,1);
	for jdx=2:m
		A(:,jdx) = A(:,jdx-1).*r;
	end
	% set up variable part of the regression matrix
	AB = zeros(nr,m,nt);
	for idx=1:nt
		for jdx=1:m
			AB(:,jdx,idx) = A(:,jdx).*B(:,idx);
		end
	end

	% solve optimisation problem
	%g_K = 1e-3;
	g_K = 2000;
	opt = optimset('TolFun',1e-15);
	LB = 0;
%	g_K = lsqnonlin(@(g_K) fopt(g_K, A, AB, beta), g_K, LB, [], opt);
	g_K = fminsearch(@(g_K) fopt(g_K, A, AB, beta), g_K, opt);
	g_K = g_K;

	% get final value
	[r x] = fopt13(g_K, A, AB, beta);

	function r = fopt(g_K, A, AB, beta)
		r = fopt13(g_K, A, AB, beta);
		r = norm(r);
		%g_K = g_K;
	end % function fopt
		
end % function prepare13

function [r x] = fopt13(g_K, A, AB, beta)
		np = size(A,2);
	nr = size(beta,1);
	nt = size(beta,2);
	x = zeros(np,nt);
	r = zeros(nr,nt);
	% solve ls problem for each time step
	for idx=1:nt
		x(:,idx) = (A - g_K*AB(:,:,idx)) \ beta(:,idx);
		% record the residual
		r(:,idx) = (A - g_K*AB(:,:,idx))*x(:,idx) - beta(:,idx);
	end
end % function fopt13


