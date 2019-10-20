% Sun Feb 15 19:04:02 CET 2015
% Karl Kastner, Berlin
%
%% Specific Discharge Method
%% upscale specific discharge to cross sectionally integrate discharge,
%%  than average
%% this method is provenly less accurate than averaging before upscaling
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
classdef SDM < HADCP_Discharge
	properties
		% coordinates
		N
		Ni
		% specific discharge weights
		w
		wi
		% velocity profile shape parameter
		alpha
%		alphai
	end
	methods
	function obj = SDM(varargin)
		obj = obj@HADCP_Discharge(varargin{:});
	%	obj.mode = mode;
	end

	% sets calibration parameter that where externally determined
	function obj = set(obj,varargin)
		if (1 == mod(length(varargin),2))
			error('input arguments must be passed in name-value pairs');
		end
		for idx=1:2:length(varargin)
			switch (varargin{idx})
			case {'N','Ni','zb','ln_z0','alpha','w','zi'}
				obj.(varargin{idx}) = varargin{idx+1};
			otherwise
				error(['unknown field',varargin{idx}]);
			end % switch
		end % for
		% interpolate values
		obj.ln_z0i  = interp1(obj.N,obj.ln_z0,obj.Ni,'linear');
		obj.z_bi    = interp1(obj.N,obj.zb,obj.Ni,'linear');
		obj.alphai  = interp1(obj.N,obj.alpha,obj.Ni,'linear');
		obj.wi      = interp1(obj.N,obj.w,obj.Ni,'linear');
	end % set

	function [c obj] = fit_(obj, h0, u_adcp, Q0);
		A = obj.regmat(h0, u_adcp);
		% estimate calibration parameter
		c = A \ Q0;
	end % fit_
	function [Q P void1 void2 obj] = predict_(obj, c, h, u_adcp)
		% specific discharges
		qi = obj.regmat(h,u_adcp);
		% scale up specific dischartes to total discharges
		Qi = bsxfun(@times,qi,1./cvec(obj.wi));
		% mean and standard error
		[Q serr] = mean_man(Qi');
%		% predict
%		Q = A*c;
		% part of constant and velocity scaled
		P = [];
		void1 = [];
		void2 = [];
	end % predict
	function [qi obj] = regmat(obj,h,u_adcp)
		z_b  = obj.zb;
		z_bi = obj.z_bi;
		ln_z0i  = obj.ln_z0i;
		alphai  = obj.alphai;
		zi      = obj.zi;
		nc      = length(h);
		ns      = length(z_b);
		ni      = length(z_bi);

		% depth N-t
		h_      = repmat(rvec(h),ns,1) + repmat(cvec(z_b),1,nc);
		hi_     = repmat(rvec(h),ni,1) + repmat(cvec(z_bi),1,nc);
		% cross sectional area (actually area/dx)
		area    = sum(h_,1)';
		% instrument height above z_b
		z       = z_bi - zi;
		z_      = repmat(z,1,nc);
		ln_z_   = repmat(log(z),1,nc);
		% roughness length
		ln_z0i_ = repmat(ln_z0i,1,nc);
		% mid-depth peak parameter
		alphai_ = repmat(alphai,1,nc);
		% velocity
		ui   = u_adcp'.*( (log(hi_) - 1 - alphai_ - ln_z0i_) ...
                                ./(ln_z_ + alphai_.*log(1-z_./hi_) - ln_z0i_));
		% specific discharge
		qi      = hi_.*ui;
%		qbar    = mean(qi,1)';
%		switch (obj.mode)
%		case {'linear'}
%			A    = [qbar];
%		case {'affine'}
%			A    = [area qbar];
%		case {'h'}
%			hbar = mean(hi_,1)';
%			A    = [qbar hbar.*qbar];
%		otherwise
%			error('');
%		end
	end % regmat
	end % methods 
end % SDM

