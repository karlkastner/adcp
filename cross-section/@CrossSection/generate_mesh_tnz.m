% Sat 13 Jan 17:45:37 CET 2018
%% generate a t-n-z mesh
function obj = generate_mesh_tnz(obj,adcp,ensmask)
	% limits
	% TODO make this a member variable
	tlim = [min(adcp.time) max(adcp.time)];
	nlim = 0.5*obj.transect.dwidth*[-1 1];
	zlim = [0 1];

	% make convex
	tlim = mean(tlim) + (1+sqrt(eps))*(tlim-mean(tlim));
	nlim = mean(nlim) + (1+sqrt(eps))*(nlim-mean(nlim));
	zlim = mean(zlim) + (1+sqrt(eps))*(zlim-mean(zlim));

	% number of points along dimensions
	switch (obj.vdimension)
	case {'tn','tnz'}
		nt   = round((tlim(2)-tlim(1))/obj.dt);
	end
	nn   = round(obj.transect.dwidth/obj.dw)-1;


	switch (class(obj.grid_n))
		case {'0d'}
			% nothing to do
		case {'SparseMesh1'}
			switch (obj.vdimension)
			case {'n','tn'}
				% nothing to do
			case {'nz'}
				obj.grid_nz = obj.grid_n;
			case {'tnz'}
				obj.grid_tnz = obj.grid_tn;
			otherwise
				error('here');
			end
		case {'Grid1'}
			switch (obj.vdimension)
			case {'n','tn'}
				% nothing to do
			case {'nz','tnz'}

				% width range
				wrange = 0.5*obj.transect.dwidth*[-1 1];
		
				% depth range
				drange = [min(obj.grid_n.val.zb), 0];
				%max   = 1.05*quantile(adcp.H(ensmask),0.95);
				%drange = [-bmax 0];
				%if (isempty(dz))
				%	dz = max(adcp.binsize);
				%end
				obj.grid_nz = Grid2('R1',wrange, 'dx1', obj.dw, 'R2', drange, 'dx2', obj.dz);
				obj.grid_nz.init();

				obj.grid_nz.fitfun = obj.grid_n.fitfun;
				obj.grid_nz.predfun = obj.grid_n.predfun;
	
				% build index		
				binmask = bsxfun(@and,adcp.mask,rvec(ensmask));

				% scale Z
				H = adcp.ens.H;
				N = adcp.bin.N(obj.cdx);	
				Z = bsxfun(@plus,adcp.bin.Z,-rvec(H));
				% mean depth at every ensemble
				% TODO apply level correction
				Nens = adcp.N(:,obj.cdx);
				d = -obj.zb(Nens);
				% scaling conserve discharge
				Z = bsxfun(@times,Z,rvec(H)./rvec(d));

				obj.grid_nz.build_index(flat(N(binmask)),flat(Z(binmask)),'i1'); %double(flat(d(binmask))));
			otherwise
				error('invalid vdimension');
			end % switch vdimension
	case {'RegularizedInterpolator1'}
	% TODO use better function
	% TODO subtract level
	% TODO use local number of cells
	hmean = -median(obj.grid_n.vali.zb);
	nz    =  round(hmean/obj.dz);

	switch (lower(obj.vdimension))
		case {'tn','n'}
			% nothing to do
		case {'nz'}
			% generate the nz-mesh
			obj.grid_nz = RegularizedInterpolator2();
			obj.grid_nz.remesh(nlim,zlim,[nn,nz]);
		case {'tnz'}
			% also generate the nz-mesh
			obj.grid_nz = RegularizedInterpolator2();
			obj.grid_nz.remesh(nlim,zlim,[nn,nz]);
	
			% tnz-mesh
			% TODO true z-mesh
			obj.grid_tnz = RegularizedInterpolator3();
			obj.grid_tnz.remesh(tlim,nlim,zlim,[nt,nn,nz]);
		otherwise
			error('invalid vdimension');
		end % switch vdimension
	otherwise
		error('invalid vmethod');	
	end % switch vmethod
end % generate_mesh_tnz

