% So 8. Nov 14:49:53 CET 2015
%% process the backscatter in 2+1D (time, across channel and along vertical)
function [obj] = process_backscatter_tnz(obj,adcp,ensmask)
	field = 'sediment_concentration';

	switch (class(obj.grid_n))
	case {'SparseMesh1'}
		switch (lower(obj.vdimension))
		case {'n','tn'}
			% nothing to do
		case {'nz'}
			error('here');
		case {'tnz'}
			S    = adcp.bin.S;
			last = adcp.ens.last;
			last = last(ensmask);
		
			% TODO check, that the ssc has been prefiltered
			val = adcp.bin.(field)(:,ensmask);

			obj.grid_tn.assignS((field),S(:,ensmask),val,last);
		otherwise
			error('here');
		end
	case {'Grid1'}
		switch (lower(obj.vdimension))
		case {'n','tn'}
			% nothing to do
		case {'nz','tnz'}
			time         = adcp.bin.time;
			val          = adcp.bin.(field);
			binmask = bsxfun(@and,adcp.mask,rvec(ensmask));

			optarg = {}; % T(:,ensmsk),N(:,ensmsk), S(:,ensmsk), ...

			% TODO scale similarly to velocity	
			obj.grid_nz.fit( ...
				      'i1' ...
				    , (field) ...
				    , val(binmask) ...
	                            , time(binmask) ...
				    , optarg{:} ...
				    );
			% invalidate out of range values
			mask_ = double(obj.grid_nz.mask);
			mask_(~mask_) = NaN;
			obj.grid_nz.val.(field) = bsxfun(@times, obj.grid_nz.val.(field),mask_);
		otherwise
			error('invalid vdimension');
		end % switch vdimension
	case {'RegularizedInterpolator1'}
		switch (lower(obj.vdimension))
		case {'n','tn'}
			% nothing to do
		case {'nz'}
			error('TODO implement');
		case {'tnz'}
				%tid    = adcp.ens.tid(:,obj.cdx);
				%msk    = logical(tid);
				binmsk = bsxfun(@and,adcp.mask,rvec(ensmask));
	
				bc = [];	
				
				% 3D is not flux averaged
				time = adcp.bin.time;
				%T = adcp.bin.T(:,cdx);
				N = adcp.bin.N(obj.cdx);
				% S is independent of cdx (transformation only in xy)
				S  = adcp.bin.S;
				val = adcp.bin.(field);
	
				lambda = [obj.lambda.t, obj.lambda.n, obj.lambda.z];
				obj.grid_tnz.init([time(binmsk) N(binmsk) S(binmsk)], ...
					[val(binmsk)], (field), lambda, [], bc);
		otherwise
				error('invalid vdimension');
		end % switch vdimension
	otherwise
		error('invalid vmethod');
	end % switch vmethod
end % process_backscatter_tnz

