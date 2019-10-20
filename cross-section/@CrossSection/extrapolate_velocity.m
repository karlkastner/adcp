% Tue Oct 28 09:18:57 CET 2014
% Karl Kastner, Berlin
%% extrapolate the velocity to the bank, bed, and surface
% TODO mark or integrate extrapolated cells
function obj = extrapolate_velocity(obj)
	switch (class(obj.grid_n))
	case {'0d'}
	case {'SparseMesh1'}
		% nothing to do
	case {'Grid1'}
	switch (lower(obj.vdimension))
	case {'nz','tnz'}

		% TODO for now, extrappolation is only implemented for U
		% TODO pass gridNZ here as an argument
		[obj.grid_nz.val.U, obj.grid_nz.eid] = obj.extrapolate_velocity_2d_STATIC(...
				  obj.grid_nz.val.U ...
				, obj.zb ...
				, obj.grid_nz.cX2 ...
				, obj.grid_nz.cXX1 ...
				, obj.grid_nz.cXX2 ...
				, obj.grid_nz.mask ...
				, []); %obj.clip );

		[obj.grid_nz.val.V, obj.grid_nz.eid] = obj.extrapolate_velocity_2d_STATIC(...
				  obj.grid_nz.val.V ...
				, obj.zb ...
				, obj.grid_nz.cX2 ...
				, obj.grid_nz.cXX1 ...
				, obj.grid_nz.cXX2 ...
				, obj.grid_nz.mask ...
				, []); %obj.clip );

		[obj.grid_nz.val.W, obj.grid_nz.eid] = obj.extrapolate_velocity_2d_STATIC(...
				  obj.grid_nz.val.W ...
				, obj.zb ...
				, obj.grid_nz.cX2 ...
				, obj.grid_nz.cXX1 ...
				, obj.grid_nz.cXX2 ...
				, obj.grid_nz.mask ...
				, []); %obj.clip );
	end
	case {'RegularizedInterpolator1'}
		% nothing to do
	otherwise
		error('unimplemented vmethod');
	end % switch
end % function extrapolate_velocity()

