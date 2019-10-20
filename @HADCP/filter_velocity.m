% 2015-03-05 14:58:41.032110278 +0100
%% filter outliers in velocity data
function obj = filter_velocity(obj,field,s)
	if (nargin()<3)
		% maximum multiple of median abslute deviations
		% TODO no magic numbers
		s = 3;
	end

	siz = size(obj.velocity.(field));
	n = prod(siz(1:2));
	for idx=1:siz(3)
		v  = obj.velocity.(field)(:,:,idx);
		v  = flat(v);
		id = (1:length(v))';
		% padd nan
		v(n+1) = NaN;
		
		% indices of 4 neighbours
		% above,below,left,right
		% (could be extended to 8 neighbours)
		nid = [id-siz(1), ...
                       id+siz(1), ...
                       id-siz(2), ...
                       id+siz(2)];
		nid(nid<=0) = n+1;
		nid(nid>n) = n+1;
		
		% value of the neighbours
		vn = v(nid);

		% remove padded nan
		v  = v(1:n);

		% robust bilinear interpolate
		vn = nanmedian(vn,2);
		
		% absolute difference
		ad = abs(v-vn);	
		mad_ = nanmedian(ad);

		% invalidate velocities that are far away from their neighbours
		fdx = ad>s*mad_;
		v(fdx) = NaN;
		obj.velocity.(field)(:,:,idx) = reshape(v,siz(1:2));
	end % for idx

%	d1     = v - [v(2,:);
 %                     0.5*(v(1:end-2,:)+v(3:end,:));
%		      v(end-1,:)];
%	%d1     = [diff(v); NaN(1,size(v,2))];
%	% TODO instead of the median a curve could be fit here
%	d2     = bsxfun(@minus,v,nanmedian(v));
%	q1     = quantile(abs(d1(:)),0.98);
%	q2     = quantile(abs(d2(:)),0.98);
%	fdx    = abs(d1) > q1 | abs(d2) > q2;
%	v(fdx) = NaN;
%	hadcp.velocity.cs(:,:,1) = v;
end % filter_velocity

