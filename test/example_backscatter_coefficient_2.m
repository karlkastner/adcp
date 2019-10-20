% Sun  8 Jul 13:38:08 CEST 2018

%f = 2.^(-1:1)*1e6;
f = [0.2,1,5]*1e6;
D = 2*linspace(0,1,1e3);

% scale to thornes representation
scale = 1/( 16/(3*pi) );


namedfigure(1,'Thorne 2012 Fig 1 b');
clf()
for idx=1:length(f)
	[ks2,fbs] = backscatter_coefficient_2(D,f(idx));
	plot(0.5*D*1e3,fbs);
	ylabel('f_{bs}')
	hold on
end
xlim([0,1e3])
ylim([0,1.4])
xlabel('Radius (um)');
legend(num2str(cvec(f),'%1.1e Hz'))

C = {'medwin-1977','thorne-2002','thorne-2008'}
f = 1.2e6;

figure(2);
clf
for idx=1:length(C)
	[ks2,fbs,chi] = backscatter_coefficient_2(D,f,C{idx});

	subplot(2,1,1)
	loglog(D,scale*ks2);
	hold on
	xlabel('D');
	ylabel('k_s^2 (thorne scale applied)');

	subplot(2,1,2)
	loglog(chi,fbs);
	hold on
	xlabel('\chi');
	ylabel('f_{bs}');
	%ylim([0,1.5]);
	ylim([10^-2,10]);
	xlim([10^-1,10^2]);
end

namedfigure(3,'Lisst Size Range');
clf();
lisst = LISST();
d_mm = 1e-3*cvec(lisst.size_um);
f    = 1.2e6;
mode = []; %'medwin-1977';
[ks2,fbs,chi,D1] = backscatter_coefficient_2(d_mm,f,mode);
D1
loglog(d_mm,ks2);
M    = 1;
[as,ass,asnu] = attenuation_coefficient(d_mm,f,M);
hold on;
loglog(d_mm,[ass,asnu,as]);
[d_mm.^0,log(d_mm)] \ log(cvec(ks2))
xlabel('D (mm)');
xlim(limits(d_mm));
legend('location','southeast','k_s^2 (no thorne scale)','\alpha_s/M')

