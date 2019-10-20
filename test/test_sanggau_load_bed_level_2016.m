% Tue 20 Dec 14:52:46 CET 2016
% Karl Kastner, Berlin

	[gridN] = sanggau_load_bed_level_2016();

if (0)
	load ../water-level/mat/kapuas-surface-level-2016.mat

	figure(1);
	clf();
if (0)
	subplot(2,1,1)
	plot(X,Y)
	hold on;
	plot([meta.xpmin meta.xpmax], [meta.ypmin,meta.ypmax]);
	axis([453188      455729    10011611    10013627]);

	plot(gridX(fdx),Y(fdx),'.');
end
	subplot(2,2,3)
	ylabel('-Depth (m)');
	plot(gridN.cX1, -gridN.val.bottom(:,1),'.-');
	xlabel('m')
end

	subplot(2,2,4)
	ylabel('Bed and surface elevation (m)');
	%plot(gridN.cX1, -gridN.val.bottom(:,1)+z0,'.-');
	plot(gridN.cX1, [gridN.val.z_b gridN.val.z_s],'.-');
	xlabel('m')

%	N = meanfilt1(cvec(N),nf);
%	D = meanfilt1(cvec(D),nf);

%	hold on
%	plot(N,-D);

