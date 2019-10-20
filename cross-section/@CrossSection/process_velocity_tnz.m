% Tue Jan 28 16:10:32 WIB 2014
% Karl Kastner, Berlin
%% process velocity data in 2+1D (time, across-section and along vertical)
function obj = process_velocity_tnz(obj, adcp, ensmask)

	switch (class(obj.grid_n))
	case {'SparseMesh1'}
		switch (lower(obj.vdimension))
		case {'n','tn'}
			% nothing to do
		case {'nz'}
			error('here');
		case {'tnz'}
			% for u,v,w
			S       = adcp.bin.S;

			U       = adcp.bin.velocity.cs(:,:,1);
			V       = adcp.bin.velocity.cs(:,:,2);
			W       = adcp.bin.velocity.cs(:,:,3);
			last    = adcp.ens.last;
			last    = last(ensmask);
		
			obj.grid_tn.assignS('Utnz',S(:,ensmask),U(:,ensmask),last);
			obj.grid_tn.assignS('Vtnz',S(:,ensmask),V(:,ensmask),last);
			obj.grid_tn.assignS('Wtnz',S(:,ensmask),W(:,ensmask),last);
		otherwise
			error('here');
		end
	case {'Grid1'}
	switch (lower(obj.vdimension))
		case {'n','tn'}
			% nothing to do
		case {'nz','tnz'}
			% fetch variables
		
			time    = adcp.bin.time;
			% N       = adcp.bin.N(obj.cdx);	
			% T       = adcp.bin.T(obj.cdx);
			U       = adcp.bin.velocity.cs(:,:,1);
			V       = adcp.bin.velocity.cs(:,:,2);
			W       = adcp.bin.velocity.cs(:,:,3);
			binmask = bsxfun(@and,adcp.mask,rvec(ensmask));
			if (0)
				Z = adcp.bin.S;
			else
				% mean depth at every ensemble
				% TODO apply level correction
				H    = adcp.ens.H;
				Nens = adcp.N(:,obj.cdx);
				d    = -obj.zb(Nens);
				% scaling conserve discharge
				U = bsxfun(@times,U,rvec(d)./rvec(H));
				V = bsxfun(@times,V,rvec(d)./rvec(H));
				W = bsxfun(@times,W,rvec(d)./rvec(H));
			end

			% average (or fit) the velocity values in the grid cells
			optarg = {};
			obj.grid_nz.fit('i1', 'Utnz', U(binmask), time(binmask), optarg{:});
			obj.grid_nz.fit('i1', 'Vtnz', V(binmask), time(binmask), optarg{:});
			obj.grid_nz.fit('i1', 'Wtnz', W(binmask), time(binmask), optarg{:});
		
			% invalidate the cells below the bottom
			val = obj.grid_nz.val;
			zb  = obj.zb;
			ZZb = repmat(zb, 1, obj.grid_nz.n2-1);
		
			obj.grid_nz.valid = obj.grid_nz.cXX2 >= ZZb;
			gmask   = (obj.grid_nz.cXX2 > ZZb);
			mask_ = double(gmask);
			mask_(~mask_) = NaN;
			val.U = bsxfun(@times,val.Utnz,mask_);
			val.V = bsxfun(@times,val.Vtnz,mask_);
			val.W = bsxfun(@times,val.Wtnz,mask_);

%			for idx=1:size(val.U,3)
%				v = val.U(:,:,idx);
%				v(~gmask) = NaN;
%				val.U(:,:,idx) = v;
%				v = val.V(:,:,idx);
%				v(~gmask) = NaN;
%				val.V(:,:,idx) = v;
%				v = val.W(:,:,idx);
%				v(~gmask) = NaN;
%				val.W(:,:,idx) = v;
%			end
		
			% write back
			obj.grid_nz.val  = val;
			%obj.grid_nz.err  = err;
			obj.grid_nz.mask = gmask;
		otherwise
			error('invalid vdimension');
		end % switch vdimension
	case {'RegularizedInterpolator1'}
		switch (lower(obj.vdimension))
		case {'n','tn'}
			% nothing to do
		case {'nz'}
			error('todo implement');
		case {'tnz'}
			% 3d processing
			% TODO no slip condition at bed
	%		bc2.Y    = 0.5*obj.transect.dwidth*[-1 1];
	%		bc2.val  = [0 0];
	%		bc2.dmax = 0.5*obj.dw;
			bc = [];
	
			% fetch bin values
	
			time = adcp.bin.time;
			%T = adcp.bin.T(:,cdx);
			N = adcp.bin.N(obj.cdx);
			% S is independent of cdx (transformation only in xy)
			S = adcp.bin.S;
			%S = S(:,obj.cdx);
			%false(adcp.ens.n,1);
			% H = adcp.bin.H;
			U = adcp.velocity.cs(:,:,1);
			V = adcp.velocity.cs(:,:,2);
			W = adcp.velocity.cs(:,:,3);
	%		UVW = adcp.velocity.cs;
	
			binmsk = bsxfun(@and,adcp.mask,rvec(ensmask));
	
			lambda = [obj.lambda.t, obj.lambda.n, obj.lambda.z];
			obj.grid_tnz.init([time(binmsk), N(binmsk), S(binmsk)], ...
				[U(binmsk), V(binmsk), W(binmsk)], 'UVW', lambda, [], bc);
			obj.grid_tnz.vali.Utnz = obj.grid_tnz.vali.UVW(:,1);
			obj.grid_tnz.vali.Vtnz = obj.grid_tnz.vali.UVW(:,2);
			obj.grid_tnz.vali.Wtnz = obj.grid_tnz.vali.UVW(:,3);
	%		obj.grid_tnz.init([time(msk),N(msk) S(msk)],adcp.ens.velocity.cs(msk,1).*H(msk),'Q',lambda,[],bc2);
	
			% standard 2d processing
	
			% fetch ens values
	
			time = adcp.ens.time;
			N = adcp.ens.N;
			N = N(:,obj.cdx);
			H = adcp.ens.H;
	%		H = H(:,obj.cdx);
	
			bc2.Y    = 0.5*obj.transect.dwidth*[-1 1];
			bc2.val  = [0, 0];
			bc2.dmax = 0.5*obj.dw;
			
			lambda   = [obj.lambda.t, obj.lambda.n];
			obj.grid_tn.init([time(ensmask), N(ensmask)], ...
					 adcp.ens.velocity.cs(ensmask,1), ...
					 'U',lambda,[],bc2);
			obj.grid_tn.init([time(ensmask), N(ensmask)], ...
					 adcp.ens.velocity.cs(ensmask,1).*H(ensmask), ...
					'Q',lambda,[],bc2);
	
			% TODO this is the determination of the cs and should go there
			bc1.X   = 0.5*obj.transect.dwidth*[-1 1];
			bc1.val = [0, 0];
			lambda  = obj.lambda.n;
			obj.grid_n.init(N(ensmask),H(ensmask),'H',lambda,[],bc1);
		otherwise
			error('invalid vdimension');
		end
	otherwise
		error('invalid vmethod');
	end %
end % process_velocity_tnz

