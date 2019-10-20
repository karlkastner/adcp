% 2014-09-02 11:13:15.009281662 +0200
% Karl Kastner, Berlin

% what is the transverse bed slope,
% if the bed is represented by a higher order polynomial?

n=9;
m = 4;
for idx=1:n
	x0 = 0.25 + 0.25*(idx-1)/(n-1);
switch(m)
case{3}
	c = [1 0 0 0;
	     1 x0 x0.^2 x0.^3;
	     0  1 2*x0 3*x0.^2;
	     1 1 1 1] \ [0 -1 0 0]';
case{4}
	c = [1 0 0 0 0;
	     1 x0 x0.^2 x0.^3 x0.^4;
	     0  1 2*x0 3*x0.^2 4*x0.^3;
	     1 1 1 1 1
	     0  1 2 3 4] ...
		 \ [0 -1 0 0 0]';
case {5}
	c = [1 0 0 0 0 0;
	     1 x0 x0.^2 x0.^3 x0.^4 x0^5;
	     0  1 2*x0 3*x0.^2 4*x0.^3 5*x0^4;
	     1 1 1 1 1 1;
	     0  1  2  3  4  5;
	     0  0  2  0  0  0] ... %0  0  2  6 12 20] ...
		 \ [0 -1 0 0 0 0]';
end
	subplot(3,3,idx);
	cla();
	x = linspace(0,1,100)';
	for jdx=1:length(c)
		A(:,jdx) = x.^(jdx-1);
	end
	y = A*c;
	plot(x,y);
	hold on
	set(gca,'xtick',linspace(0,1,11));
	c = A(:,1:2) \ y;
	plot(x,A(:,1:2)*c,'r');
	title(c(2))
end

