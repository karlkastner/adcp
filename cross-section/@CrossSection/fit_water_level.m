% Sun Feb  9 22:02:34 WIB 2014
% Karl Kastner, Berlin
%% fit water level from depth measurement
%% this works only if the ADCP is stationary
% TODO this works not well, better correct all values by their local mean depth and regress all simultaneously
function obj = fit_water_level(obj)
	% find the centre of the sampled coordinates
	xm = nanmedian(obj.bin.X);
	ym = nanmedian(obj.bin.Y);
%	xm = 0.5*(obj.xpmin + obj.xpmax);
%	ym = 0.5*(obj.ypmin + obj.ypmax);

% TODO filter zeros

	% find all points which fall within a box of r_max around it
	% and filter outliers
	% TODO no magic numbers
	r2_max = 64^2;
	R2 = (obj.ens.X - xm).^2 + (obj.ens.Y - ym).^2;
	R2 = reshape(repmat(R2,1,4),[],1);
	Z  = reshape(obj.ens.Z,[],1);
	time = reshape(repmat(obj.ens.time,1,4),[],1);
	mz = nanmean(Z);
	fdx = find(R2 < r2_max & abs(Z-mz) < 3 & isfinite(Z));
	%Z    = reshape(obj.ens.Z(fdx,:),[],1);
	Z = Z(fdx);
	time = time(fdx);

	% regress the water level in bins
	hgrid = Grid1([obj.ens.time(1) obj.ens.time(end)], obj.dt);

	if (0 == obj.hmethod)
		m = hgrid.valbinning([ones(size(Z),'single'), Z, Z.^2],time);
		% undo normalisation
		hgrid.val = hgrid.val.*m;
		% calculate the standard deviation
		n  = hgrid.val(:,1); % == m
		s  = hgrid.val(:,2);
		s2 = hgrid.val(:,3);
		% mean and standard deviation
		mean_ = s./n;
		% stadard error
		err_  = sqrt(1./(n-1).*(s2 - 1./n.*s.^2))./sqrt(n);
	else
		% regress for diurnal, semidiurnal, trend and constant
		time_ = time-time(1); % make small necessary to reduce round off errors !!!
		A = [ones(size(time_)) time_ sin(2*pi*time_) cos(2*pi*time_) sin(4*pi*time_) cos(4*pi*time_)];
		c = A \ Z;
		t = hgrid.X1(:); 
		t = t-t(1);
		A_ = [ones(size(t)) t sin(2*pi*t) cos(2*pi*t) sin(4*pi*t) cos(4*pi*t)];
		mean_ = A_*c;
		err_ = norm(A*c - Z)/sqrt(length(Z))*ones(size(mean_));
	end
	hgrid.val = [mean_ err_];
% regress also the velocity as a sine
	
% subplot(2,1,1);  plot(hgrid.X1,[mean_ mean_-err_ mean_+err_ A*c]); hold on; plot(time,Z,'.'); datetick(); subplot(2,2,3);plot(obj.ens.X,obj.ens.Y,'.','Markersize',1);  hold on; X =repmat(obj.ens.X,1,4); Y=repmat(obj.ens.Y,1,4); plot(X(fdx), Y(fdx),'r.'); axis equal; subplot(2,2,4);  hist(Z,10000)

	obj.hgrid = hgrid;
end % fit_water_level()

