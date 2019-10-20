% Mi 11. Nov 16:31:57 CET 2015
% Karl Kastner, Berlin
%% extrapolate the backscatter
% TODO mark or integrate extrapolated cells
function obj = extrapolate_backscatter(obj)
	switch (class(obj.grid_n))
	case {'Grid1'}
%{
		% fetch
		gridN = obj.gridN;
	
		% TODO, this does not extrapolate the discharge	
		% TODO extrapolate_velocity should be part of the grid class
		[gridN.val.bs] = obj.extrapolate_backscatter_1d_STATIC( ...
				gridN.cX1, gridN.val.bs, gridN.val.bottom(:,1) );
		% write back
		obj.gridN = gridN;
%}
	case {'Grid2'}
		% fetch
		gridNZ = obj.gridNZ;

		% TODO for now, extrappolation is only implemented for U
		% TODO pass gridNZ here as an argument
		for jdx=1:4
		[gridNZ.val.bs(:,:,:,jdx), obj.gridNZ.eid] = CrossSection.extrapolate_backscatter_2d_STATIC(...
				 gridNZ.val.bs(:,:,:,jdx) ...
				,obj.gridN.val.bottom(:,1) ...
				,gridNZ.cX2 ...
				,gridNZ.cXX1 ...
				,gridNZ.cXX2 ...
				,gridNZ.mask ...
				, obj.clip );
		end
%				,obj.vmethod);

		% write back
		obj.gridNZ = gridNZ;
	otherwise
		error('unimplemented vmethod');
	end % switch
end % function extrapolate_backscatter()

