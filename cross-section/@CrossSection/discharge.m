% Tue Jan 28 16:03:27 WIB 2014
% Karl Kastner, Berlin
%
%% integrate the discharge over all finite elements of the cross section
%
function [Q_t, q_tn, q_tnz, obj] = discharge(obj, t, vflag)

	if (nargin() < 2 || isempty(t))
		t = obj.t0;
	end
	if (nargin()<3)
		vflag = false;
	end

	% integrage discharge across section
	switch(class(obj.grid_n))
	case {'0d'}
		valid = obj.transect.valid;
		Q_t = interp1(obj.transect.time(valid),obj.discharge.total(valid),t,'linear');
		if (nargout() > 1)
			Q  = struct();
			Q.total = Q_t;
			Q.centre = interp1(obj.transect.time(valid),obj.discharge.centre(valid),t,'linear');
			Q.left   = interp1(obj.transect.time(valid),obj.discharge.left(valid),t,'linear');
			Q.right  = interp1(obj.transect.time(valid),obj.discharge.right(valid),t,'linear');
		end
	case {'SparseMesh1'}
		switch (obj.vdimension)
		case {'n'}
			TODO
		case {'tn','tnz'}
			N    = obj.N;
			dw   = (N(1)-N(end))/length(N);
			q_tn = var_tn(obj,'qu',t,N);
			Q_t  = dw*sum(q_tn);
		case {'nz'}
			TODO
		otherwise
			error('here');
		end
	case {'Grid1'}
		zs   = obj.level_t(t);
		zb   = obj.zb;

		% specific discharge
		[q_tn, q_tnz]  = obj.q_tn(t);

		% element width
		dw   = obj.grid_n.dx1;
		Q_t  = cvec(csdischarge(q_tn,zs,zb,dw));
		if (vflag)
			qv_tn     = obj.var_tn('qv',t);
			% TODO quick fix, improve extrapolation
			qv_tn(isnan(qv_tn))=0;
			Q_t(:,2)  = csdischarge(qv_tn,zs,zb,dw);
		end
	% TODO top, bottom, left, right
	case {'RegularizedInterpolator1'}
		switch (obj.vdimension)
		case {'n', 'nz'}
			error('TODO implement');
			%Q = obj.grid_n.mesh.integrate_1d(
			% Q = obj.grid_nZ.mesh.integrate_2d(
		case {'tn', 'tnz'}
			Q_t  = zeros(length(t),1);
			n    = obj.grid_n.mesh.np;
			q_tn = zeros(n,length(t));
			for idx=1:length(t)
				point       = [t(idx)*ones(n,1),obj.grid_n.mesh.point(:,1)];
				%point       = [t(idx)*ones(n,1),obj.grid_n.mesh.point(:,1:2)];
				% specific discharge
				q_tn(:,idx) = obj.grid_tn.mesh.interp_2d(point,obj.grid_tn.vali.qu);
			end
			% todal discharge (integral of the specific discharge)
			Q_t  = obj.grid_n.mesh.integrate_1d(q_tn);
		otherwise
			error('invalid vdimension');
		end % switch vdimension
	otherwise
		error(['invalid vmethod ' class(obj.grid_n)]);
	end % switch vmethod
end % function CrossSection::discharge()

