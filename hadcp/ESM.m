% 2015-01-25 12:53:45.489022291 +0100
% Karl Kastner, Berlin

%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
classdef ESM < HADCPM
	properties
	end
	methods
	function obj = ESM(mode)
		obj.mode = mode;
	end
	% estimate calibration parameter
	function [c obj] = fit_(obj, h0, u_adcp, q0);
		A = obj.regmat(h0,u_adcp);
		c = A \ q0;
	end
	% predict
	function [Q P obj] = predict_(obj, c, h, u_adcp)
		A = obj.regmat(h, u_adcp);
		Q = A*c;
		% part of constant and velocity scaled
		P = []; %(A(:,1)*c(1)) (A(:,2)*c(2))];
	end
	function [A obj] = regmat(obj, h, u_adcp)
		bottom = obj.bottom;
		bottomi = obj.bottomi;
		zi     = obj.zi;
		ln_z0  = obj.ln_z0;
		ln_z0i  = obj.ln_z0i;

		nc      = length(h);
		ns      = length(bottom);
		ni      = length(bottomi);

		% calculate energy slope
		z       = bottomi + zi;
		h_       = repmat(rvec(h),ni,1) + repmat(bottomi,1,nc);
		hbar = mean(h_,1)';
		sqrtS   = Constant.KAPPA./(sqrt(Constant.g*h_).*(repmat(log(z),1,nc) - repmat(ln_z0i,1,nc))).*u_adcp';
		sqrtS   = mean(sqrtS,1)';

		% calculate pseduo discharge
		h_      = repmat(rvec(h),ns,1) + repmat(bottom,1,nc);
		C_      = sqrt(Constant.g)*(log(h_) - 1 - repmat(ln_z0,1,nc))./Constant.KAPPA;
		Q_S     = sum(C_.*h_.^1.5,1)';
		
		% calibrate
		switch (obj.mode)
		case {'linear'}
			A = [Q_S.*sqrtS];
		case {'affine'}
			A = [Q_S Q_S.*sqrtS];
		case {'h'}
			A    = [Q_S.*sqrtS hbar.*Q_S.*sqrtS];
		otherwise
			error('');
		end
	end % regmat
	end % methods
end % classdef

