% Thu 13 Jul 17:28:41 CEST 2017

%Rh=double(Rh); us.m = double(us.m); A = double(calib.area); c=[Q.^0 log(Rh)]\log(Q./A); lsqnonlin(@(c) c(1)*A.*Rh.^c(2) - Q,c)
%plot(horzcat(C{:))

nf = 20;
%nf = 1;
win = hanwin(1:nf);

dw = calib.cs_(1).gridN.dx1;
%dw = 1;

H=(arrayfun(@(x) x.gridN.val.H,calib.cs_,'uniformoutput',false));
H = horzcat(H{:});
H(isnan(H)) = 0;
for idx=1:size(H,2)
	%H(:,idx) = wmeanfilt(win,H(:,idx),1)
	H(:,idx) = smooth(H(:,idx),nf);
end

A = csarea(0,-H,dw)'
Rh = csradius(0,-H,dw)'


A = double(A);
Rh=double(Rh);
Q = double(Q);
c=[Q.^0 log(Rh)]\log(Q./A); lsqnonlin(@(c) c(1)*A.*Rh.^c(2) - Q,c)

%calib.,-meanfilt1(H,nf),dw)';
%Rh = csradius(0,-meanfilt1(H,nf),dw)';

A  = csarea(rvec(calib.zs0)+meta.z_offset,wmeanfilt(win,calib.zb.A,1),dw);
Rh = csradius(rvec(calib.zs0)+meta.z_offset,wmeanfilt(win,calib.zb.A,1),dw);

A = cvec(A)
Rh = cvec(Rh)

A = double(A);
Rh=double(Rh);
Q = double(Q);
c=[Q.^0 log(Rh)]\log(Q./A); lsqnonlin(@(c) c(1)*A.*Rh.^c(2) - Q,c)


%zb = mean(wmeanfilt(win,calib.zb.A,1),2);
zb = [];
for idx=1:5
	zb(:,idx) = smooth(calib.zb.A(:,idx),nf);
end
zb = mean(zb,2);

A  = csarea(rvec(calib.zs0)+meta.z_offset,zb,dw);
Rh = csradius(rvec(calib.zs0)+meta.z_offset,zb,dw);

A = cvec(A)
Rh = cvec(Rh)

A = double(A);
Rh=double(Rh);
Q = double(Q);
c=[Q.^0 log(Rh)]\log(Q./A); lsqnonlin(@(c) c(1)*A.*Rh.^c(2) - Q,c)

plot(zb)
hline(rvec(calib.zs0)+meta.z_offset)

