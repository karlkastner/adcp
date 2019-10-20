% Wed 17 Jan 12:07:28 CET 2018
%
%% generically return value stored in field "fieldname" at time t and position N
function [val_tnz, NN, ZZ, obj] = var_tnz(obj,fieldname,t,N)
	if (nargin()<3)
		t = obj.t0;
	end
	if (nargin() < 4)
		N  = obj.N;
	end
	switch (class(obj.grid_n))
	case {'SparseMesh1'}
		switch(lower(obj.vdimension))
		case {'n','tn','nz'}
			error('not available');
		case {'tnz'}
			NN = N;
			nn = length(N);
			nt = length(t);
			N  = flat(repmat(cvec(N),1,nt));
			t  = flat(repmat(rvec(t),nn,1));
			val_tnz = obj.grid_tn.interpS(fieldname,[t,N]);
			ns      = length(obj.grid_tn.S0);
			% TODO apply surface level
			zb      = obj.var_n('zb',NN);
			ZZ       = bsxfun(@times,rvec(zb),(1-obj.grid_tn.S0));
%			ZZ      = obj.grid_tn.S0;
			%H      = obj.var_tn(H)
			val_tnz = reshape(val_tnz,ns,nn,nt);
			val_tnz = obj.extrapolate_S(NN,val_tnz);
%		sum(isnan(val_tnz(:)))
%pause
		otherwise
			error('here');
		end
	case {'0d'}
		error('unavailable');
	case {'Grid1'}
		switch (lower(obj.vdimension))
		case {'n','tn'}
			error('not available');
		case {'nz','tnz'}
			val_tnz = obj.grid_nz.predict(fieldname,t);
			if (nargin()>3)
				val_tnz = interp1(obj.grid_nz.X1,val_tnz,N,'linear');
			end
		otherwise
			error('invalid vdimension');
		end
	case {'RegularizedInterpolator1'}
		switch (lower(obj.vdimension))
		case {'n','nz','tn'}
			error('not available');
		case {'tnz'}
			%val_tnz =  NaN(obj.grid_nz.mesh.np,length(t));
			point  = [zeros(obj.grid_nz.mesh.np,1),obj.grid_nz.mesh.point];
			% TODO, interpolate in one go and reshape then
			for idx=1:length(t)
				point(:,1) = t(idx);
				[A, fdx]           = obj.grid_tnz.mesh.interpolation_matrix_3d(point);
				% TODO check for existence
				val_tnz(fdx,idx,:,:) = A*obj.grid_tnz.vali.(fieldname);
			end
		otherwise
			error('invalid vdimension');
		end
	otherwise
		error('invalid vmethod')
	end % switch
end % var_tnz

