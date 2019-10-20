% Mon 22 Jan 16:56:18 CET 2018
%% plot over time and across channel
% TODO accept N as argument
function h = plot_tn(obj,t,val,id)
	if (nargin()<4)
		id = 1;
	end
	grid = obj.grid_n;
	switch (class(grid))
	case {'SparseMesh1','RegularizedInterpolator1'}
		N = obj.N;
		N = inner2outer(cvec(N));
		t = inner2outer(cvec(t));
		[nn,tt] = meshgrid(N,t);
		quadsurf(tt,nn,val','edgecolor','none');
		set(gca,'color','none','xgrid','off','ygrid','off')
	case {'RegularizedInterpolator2'}
		grid.mesh.plot(val);
	case {'Grid1'}
		%t0 = obj.t0;
		t  = cvec(t);
		t0 = [t;2*t(end)-t(end-1)];
		N  = grid.X1;
		X  = repmat(cvec(t0),1,length(N));
		Y  = repmat(rvec(N),length(t0),1);
	
		%X = grid.XX1;
		%Y = grid.XX2;
		% top left
		X_(:,1) = flat(X(1:end-1,1:end-1));
		Y_(:,1) = flat(Y(1:end-1,1:end-1));
		%  top right
		X_(:,2) = flat(X(1:end-1,2:end));
		Y_(:,2) = flat(Y(1:end-1,2:end));
		% bottom right
		X_(:,3) = flat(X(2:end,2:end));
		Y_(:,3) = flat(Y(2:end,2:end));
		%  bottom left
		X_(:,4) = flat(X(2:end,1:end-1));
		Y_(:,4) = flat(Y(2:end,1:end-1));
		Z = flat(val');
		fdx = isfinite(Z);
		%h = patch(X_(fdx,:),Y_(fdx,:),Z(fdx));
		%h= patch(X_(fdx,:)',Y_(fdx,:)',0.*Y_(fdx,:)',Z(fdx));
		h= patch(X_(fdx,:)',Y_(fdx,:)',Z(fdx));
		%h=pcolor(X,Y,val(:,:,id).');
		h.EdgeColor='none';
		axis tight
	otherwise
		error('here');
	end
end

