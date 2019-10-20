% Thu 31 May 20:35:43 CEST 2018
% Karl Kastner, Berlin
%% return value stored in field "fieldname" at position "N" in the cross section
% TODO, optional integration over depth when z is resolved
function [val_n, N, obj] = var_n(obj,fieldname,N)
	if (nargin()<3)
		t = obj.t0;
	end
	switch (class(obj.grid_n))
	case {'0d'}
		error('unavailable');
	case {'Grid1'}
		val_n = obj.grid_n.val.(fieldname);
	case {'SparseMesh1'}
		if (nargin() < 3)
			N  = obj.N;
		end
		val_n = obj.grid_n.interp(fieldname,N);
		val_n = obj.extrapolate_n(N,cvec(val_n));
	case {'RegularizedInterpolator1'}
		val_n = obj.grid_n.vali.(fieldname);
	otherwise
		error('invalid vmethod');
	end % switch vmethod
end % var_tn

