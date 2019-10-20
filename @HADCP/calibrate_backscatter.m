% Sun Jun  8 14:07:55 WIB 2014
% Karl Kastner, Berlin
%
%% calibrate backscatter to sediment concentration by the method of Sassi
function obj = calibrate_backscatter(obj, model, M_ref, r_ref, t_ref, varargin)
	% model to convert backscatter to sediment concentration
	obj.model = model;
	obj.carg = varargin;

	% skip data gaps
	[t_ref, fdx] = select_time_slots(t_ref, obj.time);
	 % cut out valid region of reference values
	 %fdx = find(t_ref >= obj.time(1) & t_ref <= obj.time(end));
	% cut out first point if it is before first ADCP time slot (yields none in interpolation otherwise)
	if (t_ref(1) < obj.time(1))
		t_ref = t_ref(2:end);
		fdx   = fdx(2:end);
	end

	% values at calibration slots
	obj.M_ref = M_ref(fdx);
	obj.t_ref = t_ref; %(fdx);
	obj.r_ref = r_ref(fdx);

	obj.v_ref = zeros(size(obj.v,1),length(obj.t_ref),4);
	for idx=1:4
		obj.v_ref(:,:,idx)  = interp1(obj.time,obj.v(:,:,idx)',obj.t_ref,'linear')';
		%obj.v_ref  = interp1(obj.time,obj.v',obj.t_ref)';
	end
	obj.S_a_ref = zeros(size(obj.S_a,1),length(obj.t_ref),4);
	for idx=1:obj.nbeams
		obj.S_a_ref(:,:,idx)= interp1(obj.time,obj.S_a(:,:,idx)',obj.t_ref,'linear')';
	end
	obj.id_ref = interp1(obj.time,obj.idepth,obj.t_ref);

	% regress for parameters
	obj = feval(obj.regC{model}, obj, varargin{:}); %, r_k, M_ref, t_ref, varargin{:});

	obj.err    = obj.M_reg - obj.M_ref;
	obj.rmse   = sqrt(obj.err'*obj.err/(length(obj.err)-length(obj.param)));
	obj.corrMM = corr(obj.M_reg, obj.M_ref); 
	obj.codet  = 1 - obj.rmse^2/var(obj.M_ref);

	% check that error is a first order stationary process (unlikely)
	% and correct the error estimate based on a second order stationry processe
	% to obtain the corrected root mean square deviation compute f*rmse
	% to obtain the effective number of degrees of freedom compute 1/f^2*nc
	% err_i = \rho*err_{i-1} + eps
	rho = obj.err(1:end-1) \ obj.err(2:end);
	f   = sqrt( (1+rho)/(1-rho) );
	obj.autocorrErr = rho;
	obj.autocorrF   = f; 

	% parameter uncertainty
	S = obj.err'*obj.err/(length(obj.M_ref))*inv(obj.A'*obj.A);
	obj.param.std = sqrt(diag(S));
end % function calibrate()

