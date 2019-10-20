% Sun Feb 15 16:02:07 CET 2015
% Karl Kastner, Berlin
%
%% superclass for HADCP discharge estimation methods
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
classdef HADCP_Discharge < handle
	properties
		% parameter
		param = struct('val0',[]);
		% residual
	%	res0
		% rms residual
		serr0
		% coefficient of determination
		R2
		% akaike information criterion
		aic
		rat

		% bed level across section
		zb
		% roughness length across section
		ln_z0

		% spacing between z_b and ln_z0 samples
		dw = 1;


		% TODO compute these automatically

		% instrument level
		zi
		% bed level beneath HADCP bin centres
		z_bi
		% roughness length at HADCP bin centres
		ln_z0i
		% dip provile parameter at HADCP bin centres
		alphai;
		
		% method (for IVM: linear, affine)
		mode

		% parameter and error estimation by jacknifing
		jnflag = false;
	end % properties
	methods
		function obj = HADCP_Discharge(varargin)
			for idx=1:2:length(varargin)
				%switch(varargin{idx})
					obj.(varargin{idx}) = varargin{idx+1};
				%end
			end
		end % HADCPM

		function [q area sq p obj] = predict(obj, zs, u_adcp, varargin)
			if (obj.jnflag)
				[q area sq p] = obj.predict_(obj.param.val, zs, u_adcp, varargin{:}); 
			else
				[q area sq p] = obj.predict_(obj.param.val0, zs, u_adcp, varargin{:}); 
			end
		end % predict
%		function obj = goodness(obj)
%			% goodness of fit
%			res0   = obj.predict(obj.t0,obj.h0) - obj.q0;
%			serr0  = sqrt(res0'*res0/(length(obj.q0)-length(obj.param)));
%			R2     = 1 - serr0*serr0/var(obj.q0);
%			
%			obj.res0  = res0;
%			obj.serr0 = serr0;
%			obj.R2    = R2;
%		end % goodness
	end % methods
end % classdef HADCP_Discharge

