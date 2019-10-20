% 2018-01-19 17:30:56.245742137 +0800
function fit_bathymetry_2d(obj)
	% 2D bottom profile
	%
	% TODO allow for anisotropy
   	    cXX1 = obj.grid_NT.cXX1;
    	    cXX2 = obj.grid_NT.cXX2;
	    switch (upper(obj.bmethod2))
	    case {'RBM'}
            	n_order = 0;
           	 [val, err] = IRBM([N, T], H, [cXX1(:), cXX2(:)], n_order);
	    case {'IDW'}
	    	[val, err] = idw2([flat(N) flat(T)], flat(H), [cXX1(:), cXX2(:)]);
	    otherwise
		error('unimplemented bmethod2');
	    end
            obj.grid_NT.val.zb = reshape(val(:,end),obj.grid_NT.n1-1,obj.grid_NT.n2-1);
            obj.grid_NT.err.zb = reshape(err(:,end),obj.grid_NT.n1-1,obj.grid_NT.n2-1);
end

