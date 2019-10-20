% So 8. Nov 14:49:53 CET 2015
% Karl Kastner, Berlin
%
%%	process depth integrated backscatter over time t and acrross section N
%% note: backscatter is processed as flux
%%	 due to high concentration and backscatter near the bottom,
%%	the inner rpoduct of the discharge and concentration
%%	\bar u \bar c_s is not a good estimate of the 
%%       depth averaged sediment flux \bar{u c_s} 
% TODO ensmask is actually fixed for each cross section, not necessary to give as an argument
% TODO split in 2d and 3d processing
function obj = process_backscatter_tn(obj, adcp, ensmask)

	field = 'sediment_discharge';

	if (nargin() < 3)
		tid   = adcp.ens.tid(:,obj.cdx);                                
        	ensmask   = logical(tid);
	end

	switch (class(obj.grid_n))
	case {'0d'}
		error('not yet implemented');
	case {'SparseMesh1'}
		switch (lower(obj.vdimension))
		case {'n','nz'}
			% TODO u,v,w
			obj.grid_n.assign(field,adcp.ens.(field).total(ensmask));
			% obj.grid_n.assign('bs_integrated',adcp.ens.bs_integrated.total(ensmask));
			obj.grid_n.assign('rouse',cvec(adcp.sediment_concentration_profile.c(2,ensmask)));
		case {'tn','tnz'}
			obj.grid_tn.assign(field,adcp.ens.(field).total(ensmask));
			obj.grid_tn.assign('qbs',adcp.ens.qbs.total(ensmask));
			% obj.grid_tn.assign('bs_integrated',adcp.ens.bs_integrated.total(ensmask));
			obj.grid_tn.assign('rouse',cvec(adcp.sediment_concentration_profile.c(2,ensmask)));
		otherwise
			error('here');
		end
	case {'Grid1'}
		switch (lower(obj.vdimension))
		case {'n','tn','nz','tnz'}
			% build indices
			% obj.build_indices(adcp, ensmask);
	
			% dummy
			% Z    = false(adcp.ens.n,1);
			% f = @(time,T,N,Z,qbs) obj.vfunc(obj.t0,time,T,N,Z,cvec(qbs));

			% average the backscatter fluxes for each 1d-bin
			%[qbs void] = obj.grid_n.binop(...
			optarg = {}; % 
                       		     %, adcp.ens.T(ensmask) ...
                       		     %, adcp.ens.N(ensmask) ...
	                       	     %, Z(ensmask) ...
			obj.grid_n.fit('i1' ...
				     ,field ...
        	               	     , adcp.ens.(field).total(ensmask) ...
       				     , adcp.ens.time(ensmask) ...
				     , optarg{:} ...
				     );

%			obj.grid_n.fit('i1' ...
%				     ,'bs_integrated' ...
 %       	               	     , adcp.ens.bs_integrated.total(ensmask) ...
  %     				     , adcp.ens.time(ensmask) ...
%				     , optarg{:} ...
%				     );

			% write back
			%obj.grid_n.val.qbs  = qbs;

			% average rouse profile parameter for each 1d-bin
			%[rouse void] = obj.grid_n.binop( ...
			% TODO, this only works for if the profile is a Rouse class
			switch (class(adcp.sediment_concentration_profile))
			case {'Rouse_Profile'}
			obj.grid_n.fit( ...
		               'i1' ...
			     , 'rouse' ...
                             , adcp.sediment_concentration_profile.c(2,ensmask) ...
       			     , adcp.ens.time(ensmask) ...
			     , optarg{:} ...
			     );
			otherwise
			end
		otherwise
			error('invalid vmethod');
		end
		case {'RegularizedInterpolator1'}
		switch (lower(obj.vdimension))
		case {'n','nz'}
			% no shear bc
			bc.X = 0.5*obj.transect.dwidth*[-1 1];
			bc.val = [0 0];
			lambda = obj.lambda.n;
			N = adcp.ens.N;
			N = N(:,obj.cdx);

			msk = ensmask;
			obj.grid_n.init(N(ensmask),adcp.ens.(field).total(msk),field,lambda,[],bc);
			obj.grid_n.init(N(ensmask),adcp.ens.qbs.total(msk),'qbs',lambda,[],bc);
			% no bc for rouse profile, as undefined at the banks
			obj.grid_n.init(N(ensmask),adcp.rouse_profile.c(2,msk).','rouse',lambda);
		case {'tn','tnz'}
			bc.Y   = 0.5*obj.transect.dwidth*[-1 1];
			bc.val = [0, 0];
			bc.dmax = 0.5*obj.dw;
			lambda = [obj.lambda.t, obj.lambda.n];

			N = adcp.ens.N;
			N = N(:,obj.cdx);

			msk = ensmask;
			obj.grid_tn.init([adcp.ens.time(ensmask), N(msk)], ...
					adcp.ens.(field).total(ensmask),field,lambda,[],bc);
			obj.grid_tn.init([adcp.ens.time(ensmask), N(msk)], ...
					adcp.ens.qbs.total(ensmask),'qbs',lambda,[],bc);
			% no bc for rouse profile, as undefined at the banks
			obj.grid_tn.init([adcp.ens.time(ensmask) N(msk)], ...
					adcp.sediment_concentration_profile.c(2,ensmask).','rouse',lambda);
		end % switch vdimension
	otherwise
		error('invalid vmethod','process_backscatter');
	end
end % function process_backscatter

