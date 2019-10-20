% Fri Nov 14 11:56:07 CET 2014
% Karl Kastner, Berlin

function test()
m = 5;
n = 1000;

sz = linspace(0,1,m)';
H = 10;
dh = 0.25;
su = linspace(0,0.1,m)';
z0  = 1e-3;
u_s = 1e-2;
d0 = 1;
% true elevation
z = repmat(double(H - (d0:dh:0.96*H)'),1,n);
% true velocity
u = u_s/Constant.KAPPA*log(z/z0);
p = 4;

% successively increase the standard error of z
for idx=1:length(sz)
	par_(:,idx) = test_(0,sz(idx));
end
figure(1);
clf();
% TODO, into func
u_s_ = [ par_([1 3],:);
         1./par_([5 7],:)
	 par_([9 11],:)]
ln_z0_ = [ par_([2 4],:)./par_([1 3],:);
           par_([6 8],:);
           par_([10 12],:)]
subplot(2,2,1)
plot(sz,u_s_(1:p,:)');
hold on;
plot(sz,u_s_(p+1:end,:)','.');
plot(sz,ones(length(sz),1)*u_s,'k');
ylabel('u_s')
xlabel('s_z');
legend('for u','for u corrected', 'for ln z','for ln z corrected', 'non-linear for u','non linear for z');
subplot(2,2,2);
plot(sz,ln_z0_(1:p,:)'/log(2));
hold on;
plot(sz,ln_z0_(p+1:end,:)'/log(2),'.');
plot(sz,ones(length(sz),1)*log(z0)/log(2),'k');
xlabel('s_z');
ylabel('ld z_0')
legend('for u','for u corrected', 'for ln z','for ln z corrected', 'non-linear for u','non linear for z');

% successively increase the standard error of u
for idx=1:length(su)
	par_(:,idx) = test_(su(idx),0);
end
u_s_ = [ par_([1 3],:);
         1./par_([5 7],:)
	 par_([9 11],:) ]
ln_z0_ = [ par_([2 4],:)./par_([1 3],:);
           par_([6 8],:);
           par_([10 12],:)]
subplot(2,2,3)
plot(su,u_s_(1:p,:)');
hold on;
plot(su,u_s_(p+1:end,:)','.');
plot(su,ones(length(sz),1)*u_s,'k');
ylabel('u_s')
xlabel('s_u');
legend('for u','for u corrected', 'for ln z','for ln z corrected', 'non-linear for u','non linear for z');
subplot(2,2,4);
plot(su,ln_z0_(1:p,:)'/log(2));
hold on;
plot(su,ln_z0_(p+1:end,:)'/log(2),'.');
plot(su,ones(length(sz),1)*log(z0)/log(2),'k');
xlabel('s_u');
ylabel('ld z_0')
legend('for u','for u corrected', 'for ln z','for ln z corrected', 'non-linear for u','non linear for z');

pdfprint('img/simulation-z-uncertainty.eps');

%subplot(1,2,2)
%plotyy(sz,par__(1,:),sz,par__(2,:));
%hold on
%plotyy(sz,par_(1,:),sz,par_(2,:));

function par__ = test_(su,sz)
% samples of u
u_ = u + su*randn(size(u),'double');
% samples of z
z_ = z + repmat(sz*randn(1,n,'double'),size(z,1),1);
% four samples for estimating dz
dz = sz*randn(4,n,'double');
sz_ = repmat(serr(dz),size(z,1),1);
ln_z_ = log(z_);
ln_z_corr = ln_z_ - 0.5*(sz_.*sz_)./(z_.*z_);
fdx = find(z_ > 0);
fdx_ = find( z_>0 & ln_z_-0.5*sz_.*sz_./(z_.*z_) > 0);

% regress for u
A = [ln_z_(fdx), -ones(length(fdx),1)];
par = A \ u_(fdx);
par__(1:2,1) = par;

% for u corrected
A = [ln_z_corr(fdx_), -ones(length(fdx_),1)];
par = A \ u_(fdx_);
par__(3:4,1) = par;

% regression for ln z
A = [u_(fdx), +ones(length(fdx),1)];
par__(5:6,1) = A \ ln_z_(fdx);

% for ln_z corrected
A = [u_(fdx_), +ones(length(fdx_),1)];
par__(7:8,1) = A \ ln_z_corr(fdx_);

% non-linear least squares for u
ln_z = ln_z_;
par = lsqnonlin(@(par) u_(fdx) - (par(1)*ln_z(fdx) - par(1)*par(2)), [1 1]);
par__(9:10,1) = par;

% non-linear least squares for z
par = lsqnonlin(@(par) z_(fdx) - exp(1/par(1)*u_(fdx) + par(2)), [1 1]);
par__(11:12,1) = par;

end


end

