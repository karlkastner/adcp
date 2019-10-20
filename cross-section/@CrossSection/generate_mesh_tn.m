% Fri Jan 30 12:20:51 CET 2015
% Karl Kastner, Berlin
%% generate 1+1D mesh over time and across section
function obj = generate_mesh_tn(obj, adcp, ensmask)
	% limits
	tlim = [min(adcp.time) max(adcp.time)];
	nlim = 0.5*obj.transect.dwidth*[-1 1];

	% make convex
	tlim = mean(tlim) + (1+sqrt(eps))*(tlim-mean(tlim));
	nlim = mean(nlim) + (1+sqrt(eps))*(nlim-mean(nlim));

	% TODO no magic numbers
        Tlim_   = quantile(adcp.T(ensmask,obj.cdx) ,[0.02 0.98]);

	% number of points along dimensions
	nn   = round(obj.transect.dwidth/obj.dw)-1;

	% set velocity processing functions
	switch (class(obj.grid_n))
	case {'SparseMesh1'}
		tid = obj.transect.tdx;
		N   = adcp.ens.N;
		X   = adcp.ens.X;
		Y   = adcp.ens.Y;

		%obj.grid_n = SparseMesh1();
		% TODO no magic numbers
		obj.grid_n.S = 1;
		obj.grid_n.init(N(ensmask,tid));

		obj.grid_tn = SparseMesh2();
		obj.grid_tn.m = obj.grid_n.m;
		obj.grid_tn.S = [400,1];
		obj.grid_tn.init([adcp.time(ensmask),N(ensmask,tid)]);

		obj.grid_tn.assign('X',X(ensmask));
		obj.grid_tn.assign('Y',Y(ensmask));

		% TODO, this is a hack to store a mesh, should be done better
		%obj.grid_n.val.N = obj.transect.dwidth*((1:nn)'/(nn+1)-0.5);
	case {'Grid1'}
		% 1D (along cross section)
		% time is covered by functions
		obj.grid_n.R1 = nlim;
		obj.grid_n.dx1 = obj.dw;
		obj.grid_n.init();

		switch (lower(obj.vdimension))
		case {'0d'}
			error('here');
		case {'n','nz'}
			% no resolution in time, simply choose mean/median
 			obj.grid_n.fitfun  = @(v,varargin) deal(nanmedian(v),nanserr(v));
			obj.grid_n.predfun = @(c,t) c;
		case {'tn','tnz'}
			nt   = round((tlim(2)-tlim(1))/obj.dt);
			switch (obj.vmethod2)
			case {'lagrange'}
				obj.grid_n.fitfun  = @(val,t,varargin) lp_regress(obj.torder,obj.t0,t,val);
				obj.grid_n.predfun = @(t0,c,t)         lp_predict(obj.torder,obj.t0,c,t);
			case {'hermite'}
				hp2_regress_ = @(t0,t,val) hp2_regress(t0(1),t0(end),length(t0),t,val);
				obj.grid_n.fitfun = @(val,t,varargin) hp2_regress_(obj.t0,t,cvec(val));
				obj.grid_n.predfun = @(c,t) hp2_predict(obj.t0,c,cvec(t));
			case {'fourier'}
				obj.grid_n.fitfun  = @(val,t,varargin) fourier_resampled_fit(obj.T_fourier,obj.t0,t,val);
				obj.grid_n.predfun = @(c,t)   fourier_resampled_predict(obj.T_fourier,obj.t0,cvec(c),t);
				%obj.grid_n.fitfun  = @(val,t,varargin) fourier_regress(obj.T,obj.t0,t,val);
				%obj.grid_n.predfun  = @(c,t)   fourier_predict(obj.T,obj.t0,c,t);
			end
		otherwise
			error('here');
		end

		% build indices
		N1     = adcp.N(ensmask,obj.cdx);
		obj.grid_n.build_index(N1(:),'i1');
		N4     = adcp.N4(ensmask,:,obj.cdx);
		obj.grid_n.build_index(N4(:),'i4');

		% 2d-mesh for bathymetry
		% TODO no magic numbers
		% TODO also use regularizedinterpolator
      		dN     = 0.1*obj.transect.dwidth;
       		dT     = 0.25*min(dN,Tlim_(2)-Tlim_(1));
		obj.grid_nT = Grid2('R1',nlim,'dx1',dN,'R2',Tlim_,'dx2',dT);
		obj.grid_nT.init();
	case {'RegularizedInterpolator1'}
		obj.torder = [];
		obj.vfunc  = [];
		obj.pfunc  = [];
		obj.roughnessfunc = [];

		switch (lower(obj.vdimension))
		case {'n','nz'}
			%obj.grid_n = RegularizedInterpolator1();
			obj.grid_n.remesh(nlim, nn);
		case {'tn','tnz'}
			nt   = round((tlim(2)-tlim(1))/obj.dt);
			% always create depth-averaged N-mesh
			%obj.grid_n = RegularizedInterpolator1();
			obj.grid_n.remesh(nlim, nn);

			obj.grid_tn = RegularizedInterpolator2();
			obj.grid_tn.remesh(tlim,nlim,[nt, nn]);
		otherwise
			error('invalid vdimension');
		end % switch vdimension
	    otherwise
		error('invalid vmethod');
	end % switch vmethod
end % generate_mesh_tn

