% Tue 29 May 11:38:03 CEST 2018
%% plot rounds (consecutiver transects) navigated with the boat
function XYc_ = plot_rounds(obj,adcp)
	% figure();
%	scatter3(vadcp.X,vadcp.Y,d0,[],d0);
%	view(0,90);
	x0 = nanmedian(adcp.X);
	y0 = nanmedian(adcp.Y);
	t0 = adcp.time(1);
	nc = 200;
	s_ = linspace(0,1,nc);
	for idx=1:length(obj.rid)-1;
		rdx = obj.rid(idx):obj.rid(idx+1);

		
		figure(idx);
		clf();
		plot(adcp.X,adcp.Y,'.','markersize',1);
		hold on
		plot(adcp.X(rdx),adcp.Y(rdx),'.-');
		plot(adcp.X(rdx(1)),adcp.Y(rdx(1)),'og');
		axis equal

		if (0)
		figure(100+idx);
		clf
		plot([adcp.X(rdx)-x0,adcp.Y(rdx)-y0,1000*(adcp.time(rdx)-t0)],'.');
		hold on;
%		plot([adcp.X(rdx(1))-x0,adcp.Y(rdx(1))-y0,1000*(adcp.time(rdx(1))-t0)],'.');
	%	axis equal
		fdx = isfinite(adcp.X(rdx));
		s = (0:length(rdx)-1)/length(rdx);
		X_(:,idx) = interp1(s(fdx),adcp.X(rdx(fdx)),s_);
		Y_(:,idx) = interp1(s(fdx),adcp.Y(rdx(fdx)),s_);
		end
	end
	if (0)
	figure(1e3)
	subplot(2,2,1)
	plot([X_-x0]')
	subplot(2,2,2)
	plot([Y_-y0]')
%	figure(1e3+1)
	subplot(2,2,3)
	X_ = [mean(X_,2),median(X_,2)];
	Y_ = [mean(Y_,2),median(Y_,2)];
	plot([X_-x0,Y_-y0])'
	subplot(2,2,4)
	plot(adcp.X,adcp.Y,'.','markersize',1);
	hold on
	plot(X_,Y_,'.')
	end
end

