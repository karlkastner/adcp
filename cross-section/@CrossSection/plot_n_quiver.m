% Mon 22 Jan 16:55:35 CET 2018
%% plot quiver across section
function plot_n_quiver(obj,V,U)
	if (nargin()<3)
		id = 1;
	end
	grid = obj.grid_n;
	switch (class(grid))
	case {'RegularizedInterpolator3'}
		error('todo implement');
	%	grid.mesh.plot(val(:,id));
	case {'Grid1'}
		X=cvec(grid.cX1);
		%Y=grid.cXX2;
		%h=pcolor(X,Y,val(:,:,id).');
		%h.EdgeColor='none'
		%scale = 20;
			%q = quiver(X,Y,scale*V(:,:,id),scale*W(:,:,id),0);
			%q.MaxHeadSize = 0.5;
		scale = 0;
		quiver(X,zeros(size(X)),V,U,scale,'MaxHeadSize',0); %,daspect);
		%quiver_man(X,zeros(size(X)),V,U,scale); %,daspect);
		%quiver_man([X,zeros(size(X))],[V,U],scale); %,daspect);
		%quiver_man(X,zeros(size(X)),V(:,:,id),U(:,:,id),scale); %,daspect);
	otherwise
		error('here');
	end
end


