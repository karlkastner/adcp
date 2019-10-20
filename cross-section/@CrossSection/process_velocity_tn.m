% Do 30. Jul 10:12:36 CEST 2015
% Karl Kastner, Berlin
%% process the velocity data
function obj = process_velocity_tn(obj, adcp, msk)

	time = adcp.ens.time;
	T    = adcp.ens.T;
	T    = T(:,obj.cdx);
	N    = adcp.ens.N;
	N    = N(:,obj.cdx);
	H    = adcp.ens.H;
	Z    = false(adcp.ens.n,1);

	UVW  = 'UVW';
	q_C  = {'qu','qv','qw'};

	switch (class(obj.grid_n))
	case {'0d'}
		error('unavailable');
	case {'SparseMesh1'}
		switch (lower(obj.vdimension))
		case {'n','nz'}
			for idx=1:3
				obj.grid_n.assign(UVW(idx),adcp.ens.velocity.cs(msk,idx));
				obj.grid_n.assign(q_C{idx},adcp.ens.velocity.cs(msk,idx).*H(msk));

				obj.grid_n.assign('ubed',adcp.ens.velocity.bed.earth(msk,1));
				obj.grid_n.assign('vbed',adcp.ens.velocity.bed.earth(msk,2));
				obj.grid_n.assign('usurface',adcp.ens.velocity.surface.earth(msk,1));
				obj.grid_n.assign('vsurface',adcp.ens.velocity.surface.earth(msk,2));
			end
		case {'tn','tnz'}
			for idx=1:3
				obj.grid_tn.assign(UVW(idx),adcp.ens.velocity.cs(msk,idx));
				obj.grid_tn.assign(q_C{idx},adcp.ens.velocity.cs(msk,idx).*H(msk));

				obj.grid_tn.assign('ubed',adcp.ens.velocity.bed.earth(msk,1));
				obj.grid_tn.assign('vbed',adcp.ens.velocity.bed.earth(msk,2));
				obj.grid_tn.assign('usurface',adcp.ens.velocity.surface.earth(msk,1));
				obj.grid_tn.assign('vsurface',adcp.ens.velocity.surface.earth(msk,2));
			end
		otherwise
			error('here')
		end
	case {'Grid1'}
		% fetch 
		optarg = {};
		
		% velocities
		% TODO scale U by H_avg/H
		for idx=1:3
			% velocity
			obj.grid_n.fit('i1',UVW(idx),adcp.ens.velocity.cs(msk,idx),time(msk), optarg{:});
			% specific discharge
			obj.grid_n.fit('i1',q_C{idx},adcp.ens.velocity.cs(msk,1).*H(msk),time(msk),optarg{:});
		end

		obj.grid_n.fit('i1','ubed',adcp.ens.velocity.bed.earth(msk,1),time(msk),optarg{:});
		obj.grid_n.fit('i1','vbed',adcp.ens.velocity.bed.earth(msk,2),time(msk),optarg{:});
		obj.grid_n.fit('i1','usurface',adcp.ens.velocity.surface.earth(msk,1),time(msk),optarg{:});
		obj.grid_n.fit('i1','vsurface',adcp.ens.velocity.surface.earth(msk,2),time(msk),optarg{:});

		% TODO this should actually be part of the bed profile regression
		% mean depth
		%[obj.grid_n.val.H err.H] = obj.grid_n.binop('i1',f,time(msk),T(msk),N(msk),Z(msk),H(msk));
		% obj.grid_n.fit('i1','H',H(msk),time(msk),T(msk),N(msk),Z(msk));

		% TODO extrapolate_velocity should be part of the grid class
		[obj.grid_n.val.U] = obj.extrapolate_velocity_1d_STATIC( ...
				obj.grid_n.cX1, obj.grid_n.val.U, obj.level_t(obj.t0), obj.zb );
		[obj.grid_n.val.Q] = obj.extrapolate_velocity_1d_STATIC( ...
				obj.grid_n.cX1, obj.grid_n.val.qu, obj.level_t(obj.t0(1)), obj.zb );


		%[obj.grid_n.val.H] = obj.extrapolate_velocity_1d_STATIC( ...
		%		obj.grid_n.cX1, -obj.grid_n.val.H, obj.level_t(obj.t0(1)), obj.zb );
	
		% write back
		%obj.grid_n.val  = val;
		%obj.grid_n.err  = err;
	case {'RegularizedInterpolator1'}
	switch (lower(obj.vdimension))
	case {'n','nz'}
		% no shear and zero depth bc at banks
		bc.X   = 0.5*obj.transect.dwidth*[-1 1];
		bc.val = [0 0];

		lambda = obj.lambda.n;
		for idx=1:3
			obj.grid_n.init(N(msk),adcp.ens.velocity.cs(msk,idx),UVW(idx),lambda,[],bc);
			% specific discharge
			obj.grid_n.init(N(msk),adcp.ens.velocity.cs(msk,idx).*H(msk),q_C{idx},lambda,[],bc);
		end

		% TODO this should actually be part of the bed profile regression
%		obj.grid_n.init(N(msk),H(msk),'H',lambda,[],bc);
		% the regularized interpolator does not need explicit call of extrapolation
	case {'tn','tnz'}
		bc2.Y    = 0.5*obj.transect.dwidth*[-1 1];
		bc2.val  = [0 0];
		bc2.dmax = 0.5*obj.dw;
		
		lambda   = [obj.lambda.t, obj.lambda.n];
		for idx=1:3
			obj.grid_tn.init([time(msk), N(msk)],adcp.ens.velocity.cs(msk,idx),UVW(idx),lambda,[],bc2);
			obj.grid_tn.init([time(msk), N(msk)],adcp.ens.velocity.cs(msk,idx).*H(msk),q_C{idx},lambda,[],bc2);
		end

		obj.grid_tn.init([time(msk), N(msk)],adcp.ens.velocity.bed.earth(msk,1),'ubed',lambda,[],bc2);
		obj.grid_tn.init([time(msk), N(msk)],adcp.ens.velocity.bed.earth(msk,2),'vbed',lambda,[],bc2);

		obj.grid_tn.init([time(msk), N(msk)],adcp.ens.velocity.surface.earth(msk,1),'usurface',lambda,[],bc2);
		obj.grid_tn.init([time(msk), N(msk)],adcp.ens.velocity.surface.earth(msk,2),'vsurface',lambda,[],bc2);

%		bc1.X   = 0.5*obj.transect.dwidth*[-1 1];
%		bc1.val = [0 0];
%		lambda  = obj.lambda.n;
%		obj.grid_n.init(N(msk),H(msk),'H',lambda,[],bc1);
	otherwise
		error('invalid vdimension');
	end
	otherwise
		error('invalid vmethod');
	end
end % process_velocity_1d

