% Tue Mar  3 10:12:49 CET 2015
% Karl Kastner, Berlin

function discharge = test_procTrans_vele(id)

flag = 0;
%id = 2;
% note line 222 is not yet bt corrected, so 248 required
line = 250; %248;
testname = '/tmp/test1.m';
opt_ = '-dw-4-dz-1-discharge';
global vele dx dy dz

s = fileread([ROOTFOLDER,'src-external/adcptools-58/trunk/procTrans.m']);
filewrite(testname, ...
          [head(1,s), 'global vele dx dy dz;',tail(-2,head(line,s))]);

addpath('/tmp/');

meta = sanggau_metadata();
opt_ = '-dw-4-dz-1-discharge';
load([meta.filename.discharge{id}(1:end-4) opt_ '.mat']);
adcp = discharge.adcp.dat;
% keep only valid ensembles
if (flag)
	mask = false(size(adcp.timeV,1),1);
	mask(1:4525) = true;
	mask(5168:end) = true;
	adcp = ADCP.squeeze_STATIC(adcp,[],'skip',mask);
else
	mask = true(size(adcp.timeV,1),1);
end

	tid = ones(1,size(adcp.VEL,2));
	Pusr = 1e5*[   4.546551162508343   4.542045272175349
                       0.122617224468175   0.127176287124800];
	deltaZ          = 1;
	deltaN          = 4;
	DepthTransducer = 0.25;
	ShipReference = 'gps';
%	ShipReference = 'bt';
	UseExtHeading  = false;
	feval('test1',  adcp, tid, ...
			'DeltaZ', deltaZ, ...
			'DeltaN', deltaN, ...
			'DepthTransducer', DepthTransducer, ...
			'ShipReference', ShipReference, ...
			'UseExtHeading', UseExtHeading, ...
			'Pusr', Pusr);

figure(1);
clf();
figure(10);clf
l = {'U','V','W'}
for idx=1:3
	figure(1);
	subplot(2,3,idx)
	surface(vele(:,:,idx),'edgecolor','none');
	caxis([-1 1]);
	title(['Bart ' l{idx}]);
	axis ij
	subplot(2,3,idx+3)
	surface(double(discharge.adcp.velocity.earth(:,:,idx)),'edgecolor','none');
	caxis([-1 1]);
	axis ij
	title(['Karl ' l{idx}]);
	colormap jet
	figure(10);
	subplot(2,3,idx)
	v = double(discharge.adcp.velocity.earth(:,:,idx));
	d = vele(:,:,idx)-v(:,mask);
	mask_ = ~discharge.adcp.mask;
	d(~mask_(:,mask)) = NaN;
	surface(d,'edgecolor','none');
	axis ij
	subplot(2,3,idx+3)
	hist_man(d(:));
	colormap jet
end

figure(2);
clf
dx=squeeze(dx);
dy_ = squeeze(dy);
dy=1e7+squeeze(dy);
dz=-squeeze(dz);

qx = quantile(discharge.adcp.X4(:),[0.025 0.975]);
qy = quantile(discharge.adcp.Y4(:),[0.025 0.975]);
% one could only test the delta here
% this could be binned for the tilted data set and plotted / checked what yields lower variance
subplot(3,1.5,1);
plot(discharge.adcp.H4); hold on
plot(dz,'.'); hold on
subplot(3,3,3);
hist_man(discharge.adcp.H4(mask,:)-dz);
subplot(3,1.5,1+1.5);
plot(discharge.adcp.X4); hold on
plot(dx,'.'); hold on
ylim(qx);
subplot(3,3,6);
hist_man(discharge.adcp.X4(mask,:)-dx);
subplot(3,1.5,1+3);
plot(discharge.adcp.Y4); hold on
plot(dy,'.'); hold on
ylim(qy);
subplot(3,3,9);
hist_man(discharge.adcp.Y4(mask,:)-dy);
figure(3)
clf
subplot(3,1,1)
plot(discharge.adcp.H4(mask,:)-dz);
subplot(3,1,2)
plot(discharge.adcp.X4(mask,:)-dx);
subplot(3,1,3)
plot(discharge.adcp.Y4(mask,:)-dy);
figure(4)
clf();
plot(dx,dy,'.'); hold on
plot(discharge.adcp.X4,discharge.adcp.Y4);
%r = max(diff(qx),diff(qy));
%xlim(mean(qx)+[-1 1]*r);
%ylim(mean(qy)+[-1 1]*r);
xlim(qx);
ylim(qy);

figure(40);
clf();plot([nanmean(discharge.adcp.dat.nFiles.GGA.long,2)],'.'); hold on; plot(nanmean(discharge.adcp.dat.NMEAGGA.Long,2),'ro')
1.1e5*nanmean(diff([nanmean(discharge.adcp.dat.nFiles.GGA.long,2), nanmean(discharge.adcp.dat.NMEAGGA.Long,2)],[],2))
1.1e5*nanstd(diff([nanmean(discharge.adcp.dat.nFiles.GGA.long,2), nanmean(discharge.adcp.dat.NMEAGGA.Long,2)],[],2))

figure(400)
clf
X=double(discharge.adcp.X4(:));
Y=double(discharge.adcp.Y4(:));
H=discharge.adcp.H4(:);
fdx=isfinite(X.*Y.*H);
D=delaunay(X(fdx), Y(fdx));
X=X(fdx);
Y=Y(fdx);H=H(fdx);
patch(X(D)',Y(D)',H(D)','edgecolor','none');
view([0 90])
figure(401);
clf
X=dx(:);
Y=dy(:);
H=dz(:);
fdx=isfinite(X.*Y.*H);
D=delaunay(X(fdx), Y(fdx));
X=X(fdx);
Y=Y(fdx);H=H(fdx);
patch(X(D)',Y(D)',H(D)','edgecolor','none');
view([0 90])


% v=discharge.adcp.velocity.cs(:,:,1); v2 = discharge.gridNZ.val.U; clf(); hist_man(v(mask(:))); hold on; hist_man(v2(:),[],'r')
%v=lpnorm(discharge.adcp.velocity.earth(:,:,1:2),3); v2=lpnorm(discharge.adcp.velocity.cs(:,:,1:2),3); mask=discharge.adcp.mask(); clf(); hist_man([v(mask(:)) v2(mask(:))])
% 375
%clf(); scatter(dne,dse,[],'b.'); hold on; scatter(mne(1,:),mse(1,:),[],'r.'); plot([mde(1,:)' mse(1,:)' mne(1,:)']);
% 600
% subplot(2,1,1); imagesc(cvele1); caxis([-1 1]); subplot(2,1,2); imagesc(discharge.adcp.velocity.earth(:,:,1)-cvele1); caxis([-0.1 0.1])
% histogram before and after
% clf;hist_man(cvele1(fnde)); hold on; v=msh(1).vele(:,:,1); hist_man(v(:),[],'r'); v=discharge.adcp.velocity.earth(:,:,1); mask = discharge.adcp.mask; hist_man(v(mask(:)),[],'b.'); 
%
% subplot(3,1,1); imagesc(cvele1); caxis([-1 1]); subplot(3,1,2); imagesc(discharge.adcp.velocity.earth(:,:,1)-cvele1); caxis([-0.1 0.1]); subplot(6,1,5); imagesc(fnde); subplot(6,1,6); plot(isfinite(discharge.adcp.X)); ylim([-0.25 1.25])
%imagesc(lpnorm(t.msh.velocity(:,:,1:2),3) - lpnorm(t.msh.cs.velocity(:,:,1:2),3))
% 760
% imagesc(velcs(:,:,1) - 1.*discharge.adcp.velocity.cs(:,:,1)); colorbar; caxis([-1 1])
% mask=discharge.adcp.mask;
% nanmedian(cvele1(fnde(:))), nanmedian(cvele1(mask(:)))
% nanmean(cvele1(fnde(:))), nanmean(cvele1(m(:)))
% h_bt = (1-cosd(15))*mean(max(discharge.adcp.btrange,[],2)) + mean(mean(discharge.adcp.btrange,2) - min(discharge.adcp.btrange,[],2))
