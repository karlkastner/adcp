% Thu 19 Jul 09:33:34 CEST 2018
% TODO make maxi's method to be inherited
%% acoustic backscatter processing
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
classdef Backscatter < handle
	properties
		% paramter of the linear regression
		param

		% error covariance matrix of calibration
		C

		% reference distance
		R_ref  = 1;
		method = 'linear';	

		r2
		serr
		mad
		aicc
		withattenuation     = true;
		withbgconcentration = false;
		withbgattenuation   = false;
		%determine_leverage  = true;
		extended_statistics = true;
	end
	methods
		function obj =  Backscatter()
		end % constructor

		% backscatter coefficient
		function [ks2,s2_ks2] = ks2(obj)
			ks2    = 1./obj.param(1);
			s2_ks2 = obj.C(1,1)/obj.param(1)^2;
		
		end % ks2
		% normalized attenuation coefficient for sediment
		% no bias correction applied
		function [xi,s2_xi] = xi(obj,param)
			if (~obj.withattenuation)
				xi    = 0;
				s2_xi = 0;
				%warning('Attenuation was not fit and set to zero');
			else
			if (nargin() < 2)
				param = obj.param;
			end
			C = obj.C;
			xi = 5/log(10)*param(2)./param(1);	
			if (nargout()>1)
			s2_xi   =   C(1,1)/param(2)^2 ...
				  - 2*C(1,2)*param(1)/param(2)^2 ...
				  + C(2,2)*(param(1)/param(2))^2;
			end
			end
		end % xi
	
		% reparamatrized normalized attenuation coefficient
		function [b] = as(obj,param)
			if (~obj.withattenuation)
				b    = 0;
			else
				if (nargin() < 2)
					param = obj.param;
				end
				b = param(2);
			end	
		end

		function ssc0 = ssc0(obj,param)
			if (obj.withbgconcentration)
				if (nargin() < 2)
					param = obj.param;
				end
				np = 1;
				if (obj.withattenuation)
					np = np+1;
				end
				np = np+1;
				ssc0 = param(np);
			else
				ssc0 = 0;
			end
		end
		function as0 = as0(obj,param)
			if (obj.withbgattenuation)
				if (nargin() < 2)
					param = obj.param;
				end
				np = 1;
				if (obj.withattenuation)
					np = np+1;
				end
				if (obj.withbgconcentration)
					np = np+1;
				end
				np = np+1;
				as0 = param(np);
			else
				as0 = 0;
			end
		end % as0
	end % methods
end % Backscatter

