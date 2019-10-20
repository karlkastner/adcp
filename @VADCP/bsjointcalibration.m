% Wed Oct 16 11:57:47 UTC 2013
% Karl Kästner, Berlin
%% calibrate backscatter
% TODO this function silently assumes that the calibration parameters where already provided to both objects
function bsjointcalibration(obj, objC)
	I_f_1 = obj.I_f_1;
	I_f_2 = obj.I_f_2;
	beta_f_1 = obj.beta_f_1;
	beta_f_2 = obj.beta_f_2;
	M_ref = obj.M_ref;
	d_ref = obj.d_ref;
	nc = length(obj.M_ref);
	for idx=1:length(objC)
		I_f_1 = [I_f_1; objC{idx}.I_f_1];
		I_f_2 = [I_f_2; objC{idx}.I_f_2];
		beta_f_1 = [beta_f_1; objC{idx}.beta_f_1];
		beta_f_2 = [beta_f_2; objC{idx}.beta_f_2];
		M_ref = [M_ref; objC{idx}.M_ref];
		d_ref = [d_ref; objC{idx}.d_ref];
		nc = nc + length(objC{idx}.M_ref);
	end
	
	% calibrate with the joint set of calibration points
	try
	[pmean M_reg err A] = feval(obj.regC{obj.model}, ...
			I_f_1, ...
			I_f_2, ...
			beta_f_1, ...
			beta_f_2, ...
			M_ref, ...
...			[obj1.I_f_1;    obj2.I_f_1], ...
...			[obj1.I_f_2;    obj2.I_f_2], ...
...			[obj1.beta_f_1; obj2.beta_f_1], ...
...			[obj1.beta_f_2; obj2.beta_f_2], ...
...			[obj1.M_ref;    obj2.M_ref], ...
			[ones(nc,1)*obj.d_k], ...
...			[obj1.d_ref; obj2.d_ref], ...
			d_ref, ...
			[]);
	catch e
		obj.err = NaN;
		return;
	end
	% covariance matrix
	np = length(pmean);
	s2 = 1/(nc-np)*err'*err;
	pcov = s2*inv(A'*A);
	% write the joint calibration parameters into both objects
	% TODO, write also joint data back into all objects
	obj.param.mean = pmean(:)';
	obj.param.cov  = pcov; 
	obj.rmse = sqrt(err'*err/nc);
	%obj2.param.mean = pmean(:)';
	%obj2.param.cov  = pcov;
end % bsjointcalibration()

