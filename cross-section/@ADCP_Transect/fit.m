% 2018-05-27 12:09:55.385614450 +0200

-> TODO determine direction weighted with flow velocity (discharge)
-> cluster
-> centre of cross section (median coord)
-> limits of cross section

		% robust estimate of cross section
		% assume uniform distribution
		% weigh ensembles by travelled distance (boat may be moored for some time)
		w = abs(cdiff(x));
		% TODO in wquantile : also cut weights
		xq = wquantile(w,x,[p, 1-p]);
		yq = wquantile(y,x,[p, 1-p]);
		xc = mean(xq);
		yc = mean(yq);
		obj.xlim = xc + (xq-xc)/(0.5-p);
		obj.ylim = yc + (yq-yc)/(0.5-p);
		% 
		
		%X = obj.ens.X;
		%Y = obj.ens.Y;
		% TODO compute
		% exclude points, where no GPS position is available
		%fdx = find(~isnan(X.*Y));
		% method 1: regress the cross section slope and intersect
		%cc = [ones(length(X),1) X(:)] \ Y(:);
	        % direction from dominant velocity component and centre from mean ensemble coordinate
		X = obj.ens.X;
		%Y = obj.ens.Y;
		% TODO, this conflicts with processing the velocity values later
		% method 2 : select orientation orthogonal to main flow velocity
		vx = obj.vel(:,:,1);
		vy = obj.vel(:,:,2);
		dx = nanmedian(vx(obj.vdx));
		dy = nanmedian(vy(obj.vdx));
		%dx = nanmedian(reshape(obj.vel(:,:,1),[],1));
		%dy = nanmedian(reshape(obj.vel(:,:,2),[],1));
		% slope
		%cc(2) = -dx/dy; %dy/dx;
		% intersect
		%cc(1) = nanmean(Y) - nanmean(X)*cc(2);

		% TODO this is not good, compute both from n-min and nax !
		obj.xpmin = min(X);
		obj.ypmin = cc(1) + cc(2)*obj.xpmin; 
		obj.xpmax = max(X);
		obj.ypmax = cc(1) + cc(2)*obj.xpmax; 

