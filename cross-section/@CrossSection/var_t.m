% Fri 19 Jan 11:23:42 CET 2018
%% return value stored in filed "fieldname" at time t
%% cross sectionally integrated or averaged value
function [val_t, val_tn, obj] = var_t(obj,fieldname,t)
	if (nargin()<3)
		t = obj.t0;
	end
	switch (class(obj.grid_n))
	case {'0d'}
	case {'SparseMesh1'}
		[val_tn, val_tnz] = obj.var_tn(fieldname,t);
		N = obj.N;
		val_tn = obj.extrapolate_n(N,val_tn);	
		val_t = obj.dw*sum(val_tn);
	case{'Grid1'}
		[val_tn] = obj.var_tn(fieldname,t);

		% element width
		dw   = obj.grid_n.dx1;
		zs   = obj.level_t(t);
		zb   = obj.zb;

		fdx = ~isfinite(val_tn);
		if (any(any(fdx)))
			val_tn(fdx) = 0;
			warning('nan elements');
		end

		% integrate across section
		val_t  = cvec(csdischarge(val_tn,zs,zb,dw));
		%if (vflag)
		%	qv_tn     = obj.var_tn('qv',t);
		%end	
		%error('unimplemented');
	case {'RegularizedInterpolator1'}
	switch (lower(obj.vdimension))
		case {'n','nz'}
			error('not yet implemented');
		case {'tn','tnz'}
			% TODO use integration function of UnstructuredMesh
			[val_tn, val_tnz] = obj.var_tn(fieldname,t);
			val_t = sum(val_tn*obj.dw);
		otherwise
			error('here')
		end
	otherwise
		error('here')
	end
end % var_t

