% Sun Feb 15 16:27:35 CET 2015
% Karl Kastner, Berlin
%% velocity profile method
%% correct individual bin velocities for vertical velocity profile variation,
%% then averagem, then upscale to cross sectionally integrated discharge
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
classdef VPM < HADCPM
	properties
	end
	methods
	function obj = VPM(mode)
		obj.mode = mode;
	end
	function [c obj] = fit_(obj, h0, u0_adcp, Q0)
		A = obj.regmat(h0,u0_adcp);
		% calibration parameter
		c = A \ Q0;

	end
	function [Q P obj] = predict_(obj, c, h, u_adcp)
		A = obj.regmat(h,u_adcp);
		% predict
		Q = A*c;
		% part of constant and velocity scaled
		P = []; %(A(:,1)*c(1)) (A(:,2)*c(2))];
	end
	function [A obj] = regmat(obj,h,u_adcp)
		bottom = obj.bottom;
		bottomi = obj.bottomi;
		zi     = obj.zi;
		ln_z0i  = obj.ln_z0i;

		nc = length(h);
		ni = length(bottomi);
		ns = length(bottom);

		% local depth
		h_ = repmat(rvec(h),ns,1) + repmat(bottom,1,nc);
		% cs area (actually area/dx)
		area = sum(h_,1)';

		% local depth below instrument
		h = repmat(rvec(h),ni,1) + repmat(bottomi,1,nc);
		% TODO subtract the mean of h0
		hbar = mean(h,1)';
		% position of beam above bottom
		z = bottomi + zi;
		% convert velocity to depth averaged velocity
		u = u_adcp'.*( (log(h)-1-repmat(ln_z0i,1,nc))./(repmat(log(z),1,nc)-repmat(ln_z0i,1,nc)) );
		% average velocity over the profile
		ubar = mean(u,1)';
		switch (obj.mode)
		case {'linear'}
			A    = [area.*ubar];
		case {'affine'}
			A    = [area area.*ubar];
		case {'h'}
			A    = [area.*ubar hbar.*area.*ubar];
		otherwise
			error('');
		end
	end
	end % methods
end % IVM

