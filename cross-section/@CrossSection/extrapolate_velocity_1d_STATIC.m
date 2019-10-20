% Do 30. Jul 11:12:48 CEST 2015
% Karl Kastner, Berlin
%% extrapolate depth averaged velocity
% this assumes a constant step size of N
% TODO mark left, right and inner misses <- not necessary if count of good samples is available
% TODO limite range of extrapolation
function U = extrapolate_velocity_1d_STATIC(N,U,zs,zb)
	if (1 == size(U,2))
		zs = zs(1);
	end
	bottom = bsxfun(@minus,rvec(zs),cvec(zb));
	% interpolate velocity
	if (~isvector(U))
	% TODO this works only if the values in U are indeed velocities,
	% but not tidal components
		% convex hull of zeros at the bank
		n = size(U,1);
		m = size(U,2);
		U = [zeros(1,m); U; zeros(1,m)];
		row = (1:size(U,1))'*ones(1,size(U,2));
		col = ones(size(U,1),1)*(1:size(U,2));
		gdx = isfinite(U);
		fdx = ~gdx;
		interp = TriScatteredInterp(flat(row(gdx)),flat(col(gdx)),double(flat(U(gdx))),'natural');
		U(fdx) = single(interp(flat(row(fdx)),flat(col(fdx))));
		% quick fix
%		U(~isfinite(U)) = 0;
		% column wise
%		n = size(U,1);
%		N = (1:n)';
%		for idx=1:size(U,2)
%			fdx = ~isfinite(U(:,idx));
%			U(fdx,idx) = interp1([0;N(~fdx); n+1],[0;U(~fdx,idx);0],N(fdx),'linear');
%		end
		% remove hull
		U = U(2:end-1,:);
	else
		% repeat N for
	
		% left
		ldx    = find(bottom > 0,1,'first') - 1;
		% first valid velocity
		rdx    = find(isfinite(U),1,'first');
		ldx_   = rdx;
		% interpolate
		U(ldx+1:rdx-1) = interp1([ldx+0.5 rdx],[0 U(rdx)],ldx+1:rdx-1,'linear');
		% right
		rdx    = find(bottom > 0,1,'last') + 1;
		% last valid velocity
		ldx    = find(isfinite(U),1,'last');
		rdx_   = ldx;
		% interpolate
		U(ldx+1:rdx-1) = interp1([ldx rdx-0.5],[U(ldx) 0],ldx+1:rdx-1,'linear');
		% missing values
		inner  = false(size(U));
		inner(ldx_:rdx_) = true;
		gdx    =  isfinite(U) & inner;
		fdx    = ~isfinite(U) & inner;
		U(fdx) = interp1(N(gdx),U(gdx),N(fdx),'linear');
	end
end % extrapolate_velocity_1d

