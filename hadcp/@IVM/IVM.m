% Sun Feb 15 16:26:53 CET 2015
% Karl Kastner, Berlin
%% index velocity method
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
classdef IVM < HADCP_Discharge
	properties
	end
	methods
	function obj = IVM(varargin)
		obj = obj@HADCP_Discharge(varargin{:});
%		if (nargin() < 1)
%			obj.mode = 'affine';
%		else
%			obj.mode = mode;
%		end
	end % constructor
	function [c obj] = fit_(obj, zs0, u_adcp, Q0, varargin);
		A = obj.regmat(zs0, u_adcp, varargin{:});

		% estimate calibration parameter
		c = A \ Q0;
	end % fit_
	function [Q area sq P obj] = predict_(obj, c, zs, u_adcp,varargin)
		[A area sq]        = obj.regmat(zs,u_adcp,varargin{:});
		Q    = A*c;

		% should be square root or sum of squares
		sq   = abs(sq*c);

		% part of constant and velocity scaled
		P = [];
	end % predict

	% expects velocities of one ensemble in columns
	function [A area sq obj] = regmat(obj,z_s,u_adcp,htilde)
		%z_b = obj.z_b;
		nc  = length(z_s);
		%ns  = length(z_b);

		% local depth
		%hh = repmat(rvec(h),ns,1) + repmat(z_b,1,nc);
		% cs area (actually area/dx)
		%area = obj.dw*sum(hh,1)';

		area = cvec(csarea(z_s,obj.zb,obj.dw));

		% TODO subtract the mean of zs0
		hbar = mean(z_s,1)';
		ubar = nanmean(u_adcp,2);
		su   = nanserr(u_adcp,2);

		if (nargin < 4)
			htilde = z_s;
		end

		% regression matrix
		switch (obj.mode)
		case {'linear'}
			A    = [area.*ubar];
			sq   = [zeros(nc,1) area.*su];
		case {'affine'}
			A    = [area area.*ubar];
			sq   = [zeros(nc,1) area.*su];
%		case {'h'}
%			A    = [area.*ubar hbar.*area.*ubar];
		case {'h'}
			A    = [area area.*ubar area.*cvec(htilde).*ubar];
%		sq = [];
			sq = A;
		otherwise
			error('');
		end
	end % regmat
	end % methods 
end % IVM

