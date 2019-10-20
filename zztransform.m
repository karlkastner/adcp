% Thu 31 May 13:58:11 CEST 2018
%% non-linear mapping for bin coordinates when depth averages between ensembles
%% for avaraging several ensembles
%%
%%	preserve discharge w int u_avg dz = int int u dz dn = Q
%%	perserve shear stress is the same (u_avg)^2_s = mean((u_s)^2)
%%	preserve sediment transport w int u_avg c_avg dz = int int u c dz dn
%%      preserve rouse number
%%
%%	alternative : correct parameters for effects of averaging
%%
%% several approaches :
%% s-transform : z_1' = H0/H1 z_1, perserves u_bar
%%                                 does not preserve u_* (du/dz|_0)
%% clipping    : z_1' = z_1, z_1 < H0, does not preserve u_bar
%%				       unclear if H0>H1
%%                                     perserves (du/dz)_0 (u_*)
%% zz-transform : perserve both u_bar and u_
%
%% TODO this is non-monotoneous when difference in H0 and H1 is large
%
function [z0,dz0,c] = zztransform(H0,H1,z1,order)
	if (nargin() < 2)
		order = 2;
	end
	switch (order)
	case {1} % linear
		% this is also the solution for the power transform: a*z1^b = z0, a*b*0^(b-1) = 1 -> b=1
		z0 = z1*H0/H1;
	case {2}
		A  = quadraticA(H1);
       		%b  = [0; H0; 1];
       		b  = [0; 1; H0];
		c  = A \ b;
		c  = flipud(c);
		z0 = polyval(c,z1);
	case {3} % hermite
		A = hermiteA(H1);
       		b = [0; 1; H0; 1];
		c = A \ b;
		c = flipud(c);
		z0 = polyval(c,z1);
	case {4}
		% mixed
		% p  = z1/H1;
		% z0 = p.*z1*H0/H1 + (1-p).*z1;
		z0 =      0 + z1 + (H0/H1^2-1/H1).*z1.^2;
	otherwise
		error('here');
	end

function A = quadraticA(x)
	A = [ 1,0,0,0;           % f(0)
	      0,1,0,0;           % f'(x)
	      1,x,x.^2,x.^3;     % f(x)
              %0,1,2*x,3*x.^2 ... % f'(x)
	    ];
end % hermiteA

function A = hermiteA(x)
	A = [ 1,0,0,0; % f(0)
              0,1,0,0; % f'(0)
	      1,x,x.^2,x.^3; % f(x)
              0,1,2*x,3*x.^2 ... % f'(x)
	    ];
end % hermiteA

end
