% So 29. Mär 12:43:40 CEST 2015
% Karl Kastner, Berlin
%
%% string of options, for file name generation
function [optstr, obj] = optstr(obj,varargin)
        optstr = '';

	switch (class(obj.grid_n))
	case {'char'}
		optstr = [optstr,'-',obj.grid_n];
	case {'Grid1'}
		optstr = [optstr,'-',class(obj.grid_n)];
	       	optstr = [optstr,'-',lower(obj.vmethod2)];
	case {'RegularizedInterpolator1'}
		optstr = [optstr,'-',class(obj.grid_n)];
		% regularisation and averaging
		switch (lower(obj.vdimension))
		case {'n'}
			optstr = [optstr,'-lambda-n-',num2str(obj.lambda.n,'%g')];
		case {'nz'}
			optstr = [optstr,'-lambda-n-',num2str(obj.lambda.n,'%g')];
			optstr = [optstr,'-lambda-z-',num2str(obj.lambda.z,'%g')];
		case {'tn'}
			optstr = [optstr,'-lambda-t-',num2str(obj.lambda.t,'%g'), ...
						'-n-',num2str(obj.lambda.n,'%g')];
		case {'tnz'}
			optstr = [optstr,'-lambda-t-',num2str(obj.lambda.t,'%g'),  ...
						'-n-',num2str(obj.lambda.n,'%g'),  ...
						'-z-',num2str(obj.lambda.z,'%g')];
		end % vdimension
	case {'SparseMesh1'}
		optstr = [optstr,'-',class(obj.grid_n)];
		optstr = [optstr,'-m-',num2str(obj.grid_n.m)];
	otherwise
		error('here');
	end
	% TODO cat torder in case of v2 == lagrange and fourier

	switch (lower(obj.vdimension))
	case {'n'}
		optstr = [optstr,'-dn-',num2str(obj.dw)];
	case {'nz'}
		optstr = [optstr,'-dn-',num2str(obj.dw)];
		optstr = [optstr,'-dz-',num2str(obj.dz)];
	case {'tn'}
		optstr = [optstr,'-dt-',num2str(86400*obj.dt)];
		optstr = [optstr,'-dn-',num2str(obj.dw)];
	case {'tnz'}
		optstr = [optstr,'-dt-',num2str(86400*obj.dt)];
		optstr = [optstr,'-dn-',num2str(obj.dw)];
		optstr = [optstr,'-dz-',num2str(obj.dz)];
	end

	if (strcmp(lower(obj.vmethod2),'fourier'))
		optstr = [optstr,'-nfmax-',num2str(1/min(obj.T_fourier))];
	end
	%switch (class(obj.grid_n))
	%case {'Grid1'}
	%	% no reg
	%case {'RegularizedInterpolator1'}
	%case {'SparseMesh1'}
	%end % switch vmethod

        %optstr = [optstr,'-rgh-',lower(obj.velocity_profile_cls)];
        %optstr = [optstr,'-ssc-',lower(obj.sediment_profile_clf)];

	for idx=1:2:length(varargin)
        	optstr = [optstr,'-',varargin{idx}];
		val = varargin{idx+1};
		if (ischar(val))
	        	optstr = [optstr,'-',varargin{idx+1}];
		else
	        	optstr = [optstr,'-',num2str(varargin{idx+1})];
		end
	end
end % optstr

