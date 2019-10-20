% Wed 17 Jan 12:18:16 CET 2018
% Karl Kastner, Berlin
%% return values of field "fieldname" at time t and position N along cross section
%% typically depth integrated or averaged values
% TODO, optional integration over depth when z is resolved
function [val_tn, obj] = var_tn(obj,fieldname,t,N)
	if (nargin()<3)
		t = obj.t0;
	end
	if (nargin() < 4)
		N  = obj.N;
	end
	switch (class(obj.grid_n))
	case {'0d'}
		error('unavailable');
	case {'SparseMesh1'}
		switch(lower(obj.vdimension))
		case {'n','nz'}
			error('not available');
		case {'tn','tnz'}
			nn = length(N);
			nt = length(t);
			NN  = flat(repmat(cvec(N),1,nt));
			t  = flat(repmat(rvec(t),nn,1));
			val_tn = obj.grid_tn.interp(fieldname,[t,NN]);
			val_tn = reshape(val_tn,nn,nt);
			val_tn = obj.extrapolate_n(N,val_tn);
		otherwise
			error('here');
		end
	case {'Grid1'}
		switch (lower(obj.vdimension))
		case {'n','tn','tnz','nz'}
			val_tn = obj.grid_n.predict(fieldname,t);
%		case {'tnz','nz'}
%			% TODO this only works for 2d quantities
%			val_tnz = obj.var_tnz(fieldname,t);
%			val_tn = squeeze(nanmean(val_tnz,2));
		otherwise
			error('invalid vdimensions');
		end % switch vdimension
	case {'RegularizedInterpolator1'}
		switch(lower(obj.vdimension))
		case {'n','nz'}
			% time is dummy
			var_tn = obj.grid_n.vali.(fieldname);
		case {'tn','tnz'}
			N = obj.grid_n.mesh.X;
			val_tn = zeros(length(N),length(t));
			NN = repmat(cvec(N),1,length(t));
			TT = repmat(rvec(t),length(N),1);
			% TODO check for existence
			val_tn = obj.grid_tn.mesh.interp_2d([TT(:),NN(:)],obj.grid_tn.vali.(fieldname));
			val_tn = reshape(val_tn,length(N),length(t));
		otherwise
			error('invalid vdimension');
		end % switch vdimension
	otherwise
		error('invalid vmethod');
	end % switch vmethod
end % var_tn

