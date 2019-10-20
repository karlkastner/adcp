% Thu Jul 17 13:16:00 WIB 2014
% Karl Kastner, Berlin
%% fit vertical profile of the streamwise velocity
%%
%% this function will work with both ensenble data, eg. U_bin taken from ensembles,
%% as well as gridded data, (U_bin taken from the velocity grid)
%%
%% input
%% cs    : struct        : cross section averaged data
%% U_bin : [nrow x ncol] : vertical profiles of stream wise velocity
%% Z_bin : [nrow x ncol] : positions of bin above bottom for each element in U_bin
%% ens.N : [ncolx1]      : position of each column of U in along the cross section
%% ens.H : [ncolx1]      : depth of each column of U
%% ens.sH  : [ncolx1]     : std of depth at each colum of U
%% ens.ldx : [ncolx1]    : last valid sample in column of U
%% dw_z0   : scalar    : grid cell size for grid_n
%% obj.roughnessmethod           : method to use for the computation
%% output:
%% grid_n : struct       : function of u_s and z_0 along cross section
%% us_ens, ln_z0_ens, U_ens : local estimates for input ensembles/grid columns
%%            not returned by every obj.roughnessmethod

% TODO cs or sw direction ?
function [obj] = fit_vertical_profile_of_velocity(obj, adcp, ensmask)


    %profile = adcp.velocity_profile;
    %profile_cls = class(adcp.velocity_profile);
    param   = adcp.velocity_profile.param;

   % remove outliers
   % TODO, two pass / weighted regression to remove outliers
   if (0)
	while (1)
		p = param(1,:);
		% TODO make this a function
		fdx = find(isfinite(p));
		p = p(fdx);
		d = abs([NaN, p(1:end-1)]-p);
		d = max(d,abs([p(2:end),   NaN]-p));
		s = quantile(d,normcdf(1/2));
		m = maxnnormals(0,s,length(d));
		fdx_invalid = fdx(abs(d)>m);
		if(isempty(fdx_invalid))
			break;
		end
		param(:,fdx_invalid) = NaN;
	end

if (0)
	figure();
	clf
	subplot(2,2,1)
	plot(p0(1,:))
	hold on
	plot(param(1,:))
	subplot(2,2,2)
	plot(adcp.ens.N,p0,'.');
	hold on
	plot(adcp.ens.N,param(1,:),'.');
	subplot(2,2,3)
	scatter3(adcp.N,adcp.time,param(1,:)',[],param(1,:),'.');
	drawnow
end
   end	% outlier removal


   msk = logical(adcp.ens.tid(:,obj.cdx));

   switch (class(obj.grid_n))
   case {'SparseMesh1'}
	name   = {'a','b','c','d'};
	switch (lower(obj.vdimension))
	case {'n','nz'}
		for idx=1:size(adcp.velocity_profile.param,1)
			% function for convenience in sparsemesh
			obj.grid_n.assign(name{idx},adcp.velocity_profile.param(idx,ensmask)');
			% TODO regularise here
        		cs_param(idx,:) = obj.grid_n.val.(name{idx});
		end
	case {'tn','tnz'}
		for idx=1:size(adcp.velocity_profile.param,1)
			obj.grid_tn.assign(name{idx},adcp.velocity_profile.param(idx,ensmask)');
        		cs_param(idx,:) = obj.grid_tn.val.(name{idx});
		end
	end % switch name
    % recover the physical parameters
    % only works if points are lagrangian/hermite polynomials entirely derived from point values
    %siz       = size(cs_param);
    %cs_param_ = reshape(cs_param,siz(1)*siz(2),[]).';
%    profile = feval(profile_cls,cs_param);
    profile = objcopy(adcp.velocity_profile);
    profile.S2 = [];
    profile.serr = [];
    profile.param = cs_param;

	siz = [1,length(cs_param)];
   case {'Grid1'}
	    % count ensembles per bin (includes invalid ensembles)
	    % TODO, this should go where the index is build
	    % [obj.grid_n.val.n] = cvec(obj.grid_n.binop('i1', @sum, rvec(ones(size(time)))));

	    % dummys
   	    time  = adcp.ens.time;

	% pre-filter
   for idx=1:size(param,1)
	nf = 11;
	fdx = isfinite(param(idx,:));
	param(idx,fdx) = medfilt1(double(param(idx,fdx)),nf);
   end



	    % bin the profile parameters
	    optarg = {};
	    for idx=1:adcp.velocity_profile.np
			[cs_param(:,:,idx), void] = obj.grid_n.binop('i1', ...
						       obj.grid_n.fitfun, ...
						       cvec(param(idx,ensmask)), ...
						       time(ensmask) );
	    end
	    % smooth shear velocity
	    % this has to precede computation of roughness length
	    if (obj.dwf > 0)
		    	nf  = round(obj.dwf/obj.grid_n.dx1);
			win = kaiserwin(-nf:nf);
			[cs_param(:,:,1) s_] = wmedfilt(win,cvec(double(cs_param(:,:,1))));
	    end
    % recover the physical parameters
    % only works if points are lagrangian/hermite polynomials entirely derived from point values
    siz       = size(cs_param);
    cs_param_ = reshape(cs_param,siz(1)*siz(2),[]).';
%    profile = feval(profile_cls,cs_param_);
    profile = objcopy(adcp.velocity_profile);
    profile.S2 = [];
    profile.serr = [];
    profile.param = cs_param_;
   case ('RegularizedInterpolator1')
    switch (lower(obj.vdimension))
    case {'n','nz'}
	name   = {'a','b','c','d'};
	lambda = [obj.lambda.us_scale 1 1]*obj.lambda.n;
	for idx=1:adcp.velocity_profile.np
		if (1==idx)
			% no slip bc for shear velocity
			bc.X   = 0.5*obj.transect.dwidth*[-1 1];
			bc.val = [0 0];
		else
			bc = [];
		end
		obj.grid_n.init(adcp.ens.N(msk,obj.cdx),cvec(param(idx,msk)), name{idx},lambda(idx),[],bc);
		cs_param(:,:,idx) = obj.grid_n.vali.(name{idx});
	end
    % recover the physical parameters
    % only works if points are lagrangian/hermite polynomials entirely derived from point values
    siz       = size(cs_param);
    cs_param_ = reshape(cs_param,siz(1)*siz(2),[]).';
%    profile = feval(profile_cls,cs_param_);
    profile = objcopy(adcp.velocity_profile);
    profile.S2 = []; profile.serr = [];
    profile.param = cs_param_;
    case {'tn','tnz'}
	name   = {'a','b','c','d'};
	lambda = [obj.lambda.us_scale 1 1];
	l      = [obj.lambda.t obj.lambda.n];
	% for each parameter
	for idx=1:adcp.velocity_profile.np
		if (1==idx)
			bc.Y    = 0.5*obj.transect.dwidth*[-1 1];
			bc.val  = [0, 0];
			bc.dmax = 0.5*obj.dw;
		else
			bc = [];
		end
		N = adcp.ens.N;
		obj.grid_tn.init([adcp.ens.time(msk),N(msk,obj.cdx)], ...
			     cvec(param(idx,msk)),name{idx},lambda(idx)*l,[],bc);
		cs_param(:,:,idx) = obj.grid_tn.vali.(name{idx});
	end
    % recover the physical parameters
    % only works if points are lagrangian/hermite polynomials entirely derived from point values
    siz       = size(cs_param);
    cs_param_ = reshape(cs_param,siz(1)*siz(2),[]).';
    %profile = feval(profile_cls,cs_param_);
    profile = objcopy(adcp.velocity_profile);
    profile.S2 = []; profile.serr = [];
    profile.param = cs_param_;
    otherwise
	error('invalid vdimension');
    end
    otherwise
	error('invalid vmethod');
    end

    [v, s]   = profile.shear_velocity();
    v       = reshape(v,siz(1:2));
    val.u_s = v;
    err.u_s = s;

    [v, s, b]  = profile.ln_z0();
    v          = reshape(v,siz(1:2));
    val.ln_z0  = v;
    err.ln_z0  = s;
    bias.ln_z0 = b;

    switch (class(adcp.velocity_profile))
    case {'Log_profile'}
	% nothing to do
	val.perturbation = [];
	err.perturbation = [];
    case {'Log_profile_with_dip'}
   	[v, s] = profile.dip_parameter();
        v     = reshape(v,siz(1:2));
	val.perturbation = v;
	err.perturbation = s;
    case {'Log_profile_with_wake'}
   	[v, s] = profile.perturbation_parameter();
        v      = reshape(v,siz(1:2));
	val.perturbation = v;
	err.perturbation = s;
    otherwise
	% nothing to do
	val.perturbation = [];
	err.perturbation = [];
    end

   switch (class(obj.grid_n))
   case {'SparseMesh1'}
%        p = feval(class(adcp.velocity_profile));
%	name = {'shear_velocity','roughness_length','perturbation_parameter'}
	%for idx=1:length(name)
	%	obj.grid_n.val.(name{idx}) = p.(name);
	%end
   	   switch (lower(obj.vdimension))
	   case {'n','nz'}
		obj.grid_n.val.u_s          = val.u_s;
		obj.grid_n.val.ln_z0        = val.ln_z0;
		obj.grid_n.val.perturbation = val.perturbation;
	   case {'tn','tnz'}
		obj.grid_tn.val.u_s          = val.u_s;
		obj.grid_tn.val.ln_z0        = val.ln_z0;
		obj.grid_tn.val.perturbation = val.perturbation;
	   end
   case {'Grid1'}
	obj.grid_n.val.u_s    = val.u_s;
	obj.grid_n.val.ln_z0  = val.ln_z0;
	obj.grid_n.val.perturbation = val.perturbation;
	obj.grid_n.err.u_s    = err.u_s;
	obj.grid_n.err.ln_z0  = err.ln_z0;
	obj.grid_n.val.bias.ln_z0 = bias.ln_z0;
	obj.grid_n.err.perturbation = err.perturbation;
   case {'RegularizedInterpolator1'}
   	   switch (lower(obj.vdimension))
	   case {'n','nz'}
		obj.grid_n.vali.u_s   = val.u_s;
		obj.grid_n.vali.ln_z0 = val.ln_z0;
		obj.grid_n.vali.perturbation = val.perturbation;
	   case {'tn','tnz'}
		obj.grid_tn.vali.u_s   = val.u_s;
		obj.grid_tn.vali.ln_z0 = val.ln_z0;
		obj.grid_tn.vali.perturbation = val.perturbation;
	   otherwise
		error('invalid vdimension')
	   end
   otherwise
	error('invalid vmethod')
   end

end % fit_vertical_profile

