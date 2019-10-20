% Sun Feb 15 16:02:07 CET 2015
% Karl Kastner, Berlin
%% fit the model parameter for HADCP discharge prediction,
%% estimate errors with the Jacknife method
function obj = fit(obj, zs0, u0_adcp, q0, varargin)
		%, zi, ...
		%z_b, z_bi, ln_z0, ln_z0i, alphai, varargin)

%	obj.z_b     = z_b;
%	obj.z_bi    = z_bi;
%	obj.zi      = zi;
%	obj.ln_z0   = ln_z0;
%	obj.ln_z0i  = ln_z0i;

%			param.estimate(h0,u0_adcp,q0);

	% calibrate parameters with Jacknife
	if (obj.jnflag)
		obj.param = Jackknife(@obj.fit_);
	else
		obj.param = struct();
	end
	% calibrate parameters by least squares reqression
	obj.param.val0 = obj.fit_(zs0,u0_adcp,q0,varargin{:});

	% standard error of calibration
	%f = @obj.predict_;
	%[val bias serr0] = param.apply(f,h0,u0_adcp);
	%[val bias serr0] = param.apply(@obj.predict,h0,u0_adcp);
	val       = obj.predict(zs0,u0_adcp,varargin{:});
	res0      = val - q0;
	serr0     = sqrt(res0'*res0/(length(q0)-length(obj.param.val0)));
	R2        = 1 - serr0*serr0/var(q0);
	obj.aic	  = akaike_information_criterion(serr0,length(q0),length(obj.param.val0));

%	[Q] = obj.predict(zs0,u0_adcp,varargin{:});

	% write back
	% affine part
%	obj.rat   = NaN;
%	obj.res0  = res0;
	obj.serr0 = serr0;
	obj.R2    = R2;

	% estimate by linearisation
	% obj.goodness();
end

