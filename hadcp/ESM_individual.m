% 2015-01-25 19:53:47.677034948 +0800
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
classdef ESM_individual < handle
	properties
		affine;
		c
	end
	methods
	function obj = ESM_individual(affine)
		obj.affine = affine;
	end

	% estimate calibration parameter
	function obj = calibrate(obj,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func,bottom_,ln_z0_)
		nc = length(Q0);
		ns = size(u0,1);
		% calculate friction slope
		z = bottom + zi;
		h = repmat(rvec(l0),ns,1) + repmat(bottom,1,nc);
		sqrtS = Constant.KAPPA*u0./(sqrt(Constant.g*h).*(repmat(log(z),1,nc) - repmat(ln_z0,1,nc)));
		% calculate pseduo discharge
		ns_     = length(bottom_);
		h_      = repmat(rvec(l0),ns_,1) + repmat(bottom_,1,nc);
		C_      = sqrt(Constant.g)/Constant.KAPPA*(log(h_) - 1 - repmat(ln_z0_,1,nc));
		Q_S     = sum(C_.*h_.^1.5)';
		% calibrate for each point in the cross section an individual constant
		for idx=1:ns
			if (obj.affine)
				A = [Q_S Q_S.*sqrtS(idx,:)'];
			else
				A = Q_S.*sqrtS(idx,:);
			end
			obj.c(:,idx) = A \ Q0;
		end
	end
	% predict
	function [Qp Q_i obj] = predict(obj,zi,l0,Q0,area,q0,u0,ln_z0,bottom,func,bottom_,ln_z0_)
		nc = size(u0,2);
		ns = size(u0,1);
		% friction slope over measured range
		z = bottom + zi;
		h = repmat(rvec(l0),ns,1) + repmat(bottom,1,nc);
		sqrtS = Constant.KAPPA*u0./(sqrt(Constant.g*h).*(repmat(log(z),1,nc) - repmat(ln_z0,1,nc)));

		% Q/S over entire cross section
		nc_     = length(l0);
		ns_     = length(bottom_);
		h_      = repmat(rvec(l0),ns_,1) + repmat(bottom_,1,nc_);
		C_      = sqrt(Constant.g)/Constant.KAPPA*(log(h_) - 1 - repmat(ln_z0_,1,nc_));
		Q_S     = sum(C_.*h_.^1.5)';
		% predict, i.e. scale by friction slope
		% for all measurement points individually
		Q_i = zeros(size(sqrtS));
		for idx=1:ns
			Q_i(idx,:) = (Q_S.*([ones(nc,1) sqrtS(idx,:)']*obj.c(:,idx)))';
		end
		% TODO use special averaging function
		Qp = mean(Q_i)';
	end
	end % methods
end % classdef

