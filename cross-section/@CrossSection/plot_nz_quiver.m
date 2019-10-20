% Mon 22 Jan 16:55:35 CET 2018
%% quiver plot of velocity across section
function plot_nz_quiver(obj,V,W,id,scale,sig)
	if (nargin()<4)
		id = 1;
	end
	if (nargin() < 5 || isempty(scale))
		scale = 10;
	end
	if (nargin() < 6)
		sig = 1;
	end
	grid = obj.grid_nz;
	switch (class(grid))
	case {'double'} % TODO quick fix 'sparsemesh'}
		% TODO make a function
		% TODO apply surface level
		zb      = obj.var_n('zb',obj.N);
		N       = sig*obj.N;
		Z       = bsxfun(@times,rvec(zb),(1-obj.grid_tn.S0));
		N       = ones(size(Z,1),1)*N;
		% note : v was already been turned, no need to turn again
		quiver_man(N,Z,V(:,:,id),W(:,:,id),scale,daspect);
	case {'RegularizedInterpolator3'}
		error('here');
		%grid.mesh.plot(val(:,id));
	case {'Grid2'}
		X=grid.cXX1;
		Y=grid.cXX2;
		%h=pcolor(X,Y,val(:,:,id).');
		%h.EdgeColor='none'
			%q = quiver(X,Y,scale*V(:,:,id),scale*W(:,:,id),0);
			%q.MaxHeadSize = 0.5;
		quiver_man(X,Y,V(:,:,id),W(:,:,id),scale,daspect);
	otherwise
		error('here');
	end
end


