% Mon 22 Jan 16:56:18 CET 2018
%% plot along n and z
function h = plot_nz(obj,N,val,id,p_resample)
	if (isempty(N))
		N = obj.N;
	end
	if (nargin()<4 || isempty(id) )
		id = 1;
	end
	if (nargin()< 5 || isempty(p_resample))
		p_resample = 5;
	end
	grid = obj.grid_nz;
	switch (class(grid))
	case {'double'} % 'sparsemesh'} TODO quick fix
		% TODO apply surface level
		N_ = obj.N;
		if (sign(N(1)) ~= sign(N_(1)))
			zb      = obj.var_n('zb',-N);
		else
			zb      = obj.var_n('zb',N);
		end
		S  = obj.grid_tn.S0;
		Si = linspace(0,1,round(length(S)*p_resample))'; 
		Si = Si*ones(1,round(length(N)*p_resample)); %ones(size(z,1),1)*Ni;
		Ni = innerspace(N(1),N(end),round(length(N)*p_resample)); 
		Ni = ones(size(Si,1),1)*Ni;
		S_ = S*ones(1,length(N));
		N_ = ones(length(S),1)*rvec(N);
		% TODO use makikma in newer versions
		vi = interpn(S,N,val(:,:,id),Si,Ni,'spline');

		% centres
		%zi  = bsxfun(@times,rvec(zb),(1-Si));
		zbi  = interp1(N(1,:),rvec(zb),Ni(1,:),'linear');
		zi  = (1-Si(:,1))*zbi;
		zi  = inner2outer(inner2outer(zi)')';
		Ni  = inner2outer(inner2outer(Ni)')';
		%v  = inner2outer(inner2outer(ones(size(val(:,:,id))))')';
		quadsurf(Ni,zi,vi,'edgecolor','none');
		view(0,90)
%pause
%		quadsurf(N,z,val(:,:,id),'edgecolor','none');
		%[id1,id2] = meshgrid(1:size(N,1),1:size(N,2));
		%n = size(N);
		% corners
	%	patch([flat(N(1:end
%		surface(ones(size(z,1),1)*N,z,val(:,:,id),'edgecolor','none');
	case {'RegularizedInterpolator3','RegularizedInterpolator2'}
		grid.mesh.plot(val(:,id));
	case {'Grid2'}
		obj.grid_nz.plot(val(:,:,id));
	otherwise
		error('here');
	end
end % plot_nz

