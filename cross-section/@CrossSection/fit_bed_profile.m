% Wed Oct 29 09:52:33 CET 2014
% Sun Jan  5 14:40:34 WIB 2014
% Karl Kastner, Berlin
%
%% fit the bed profile, has to precede n-z meshing of the cross-section
function [obj] = fit_bed_profile(obj,adcp,msk)
	% depth below mean level
	zb1 = adcp.ens.zb(msk);
	zb4 = adcp.ens.zb4(msk);
%	zb1 = adcp.ens.level(msk)  - adcp.ens.H(msk);
%	zb  = adcp.ens.level4(msk) - adcp.ens.H4(msk,:);
	N4   = adcp.N4(msk,:,obj.cdx);
	T    = adcp.T4(msk,:,obj.cdx);
	N1   = adcp.N(msk,obj.cdx);
	
	% bin the depth values into the grid cells
	switch (class(obj.grid_n))
	case {'SparseMesh1'}
		obj.grid_n.assign('zb',zb1);
	case {'Grid1'}
		%obj.grid_n.fit('i4','zb',zb);
		optarg = {};
		[obj.grid_n.val.zb, obj.grid_n.err.zb] = obj.grid_n.binop('i4', obj.bfunc, zb4, optarg{:});
		obj.extrapolate_bed_profile();
	
		if (0)	
			% IDW
			func = @(x,y) idw1(x,y,0);
			% determine bottom profile by means of radial basis functions
    			cX1 = grid_n.cX1;
    			cY1 = zeros(size(cX1));
			r0 = 0.1*width;
    			[val, err, R2]   = IRBM([N, T], H, [cX1(:) cY1(:)], 0, r0);
		end
	case {'RegularizedInterpolator1'}
		if (~isempty(obj.topofbank))
			bc.X   = 0.5*obj.dwidth*[-1 1];
			bc.val = obj.topofbank*[1 1];
		else
			bc = [];
		end
		obj.grid_n.init(N1(:),zb1(:),'zb',obj.lambda.n,[],bc);
	otherwise
		error('here');
	end
end % fit_bed_profile

