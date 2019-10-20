% Sun Feb 15 16:26:53 CET 2015
% Karl Kastner, Berlin
%% Index velocity method of Horizontal ADCP data
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
classdef HADCP_IVM < HADCPM
	properties
	end
	methods
	function obj = HADCP_IVM()
	end
	function [c obj] = fit_(obj, h0, u_adcp, Q0);
		A = obj.regmat(h0, u_adcp);
		% estimate calibration parameter
		c = A \ Q0;
	end % fit_
	function [Q P obj] = predict_(obj, c, h, u_adcp)
		A    = obj.regmat(h,u_adcp);
		Q    = A*c;
		% part of constant and velocity scaled
		P = [(A(:,1)*c(1)) (A(:,2)*c(2))];
	end % predict
	function [A obj] = regmat(obj,h,u_adcp)
		bottom = obj.bottom;
		nc = length(h);
		ns = length(bottom);

		% local depth
		h = repmat(rvec(h),ns,1) + repmat(bottom,1,nc);
	
		% cs area (actually area/dx)
		area = sum(h,1)';
		hbar = mean(h,1)';
		ubar = mean(u_adcp,2);
		% regression matrix
		A    = [area.*ubar hbar.*area.*ubar];
	end % regmat
	end % methods 
end % IVM

