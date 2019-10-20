% Do 5. Nov 12:30:41 CET 2015
% Karl Kastner, Berlin
% Jan 25 19:53:13 WIB 2014
% Karl Kastner, Berlin
%
%% Level-3 ADCP data processing, projection to cross section and integration/averaging
%
% 0	Test: flip left and right end of cross section and check for identical processing (especially sign change of u and v)
% 1	determine bottom flow, side flow for regtn
% 2	prediction function : fix zero in last column for constant model
% 3	discard outliers in depth, outliers in vel and bs ?
% 4	recheck 0d-processing
% 5	verify: difference between bottom track and gps
% 3	XY-interpolation, only interpolate, when distance is not larger than threshold, revert to bottom track
% 10	true location on bottom (correct x and y by heading and theta for each beam)
% 11	todo automated transect detection
% 12	transform to s-coords before binning, reason : if depth of two ensembles differs, than velocity is smeared at the bottom without interpolating to s-coords
% 100	exclude samples, where boat does not go parallel to cross section (45 deg level)


%	- implement fast prediction
%	- implement hermite elements with support points (no derivatives) only

% TODO distinguish between gridN and gridNR with isa and do not have separate member variables
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
classdef CrossSection < handle
	properties % (SetAccess = protected)
		% time steps for regression and integration
		t0;

		% cross section grid 1D (N direction)
		grid_n

		% cross section grid 2d (N,Z)
		grid_nz

		% cross section grid 1+1D (T,N)
		% only used for vmethod=regularizedinterpolator,
		% for vmethod=grid, time is resolved by function and t-grids are not used
		grid_tn

		% cross section grid 2+1d (T,N,Z)
		grid_tnz

		% 2D grid in NT direction (top view, for bathymetry only)
		grid_nT


		% cut off radius of bottom length grid
		R_b  = [];

		%
		% channel geometry
		%

		% cut off distance in streamwise direction
		T_max = [];
		T_rel = [];

		% cross section index
		cdx = 1;

		% TODO make a function
		%qbs;

		% water level input
		level = struct('time',[],'val',[]);

		% estimation fuctions (can later be replaced by harmonic analysis)
		% median does not lead to significant differences,
		% however, mean is required to make jacknife error estimation work
		% order 0 == mean
		%vfunc = @tmedian;
		%vfunc = @(t,x) deal(nanmedian(x),nanserr(x));
		% TODO should also accept T (streamwise coordinate, this allows for idw, polynomial or kriging/covriance models)
		bfunc = @(x) deal(nanmean(x),nanserr(x));

		% method for velocity processing
		% 0d     : traditionall processing, no mesh
		% Grid1D : depth averaged, multiple time steps supported by function
		% Grid2D : depth resolved, multiple time steps supported by function
		% regnD  : depth averaged, time ageraged
		% regtnD  : depth averaged, time discretised
		%vmethod    = []; % 'Grid1';
		vdimension = 'n';
		vmethod2   = 'hermite'; % TODO merge with vfunc
		torder     = [];

		% TODO rename into vefunc and vpfunct
		vfunc = [];
		pfunc = [];

		%bclass   = 'RegularizedInterpolator1';
		bmethod   = 'FDM';
		bmethod2d = 'IDW';

		% roughness length estimation data sources
		% ENSEMBLE : ensemble data
		% VGRID    : gridded data
		roughnesssource = 'ENSEMBLE';
		
		% roughness lengths estimation methods
		% OLS_U : ordinary least squares for velocity (u = f(z) + eps)
		% OLS_Z : ordinary least squares for depth (z = f(u) + eps)
		% Z_U   : depth averaged method
		% LOG
		% DIP
		%velocity_profile_cls = 'Log_profile_with_wake';
		%sediment_profile_cls = 'Rouse_Profile';
		%roughnessfunc = @median_man 
		%roughnessfunc = @(t,x) tmedian(t0,t,x), ...
		roughnessfunc =  @(t0, time, T, N, Z, val) median_man(val);
		% {@(t,x) median_man(x), ...
                %                @(t,x) sum(x) };

		% thikonov regularisation parameter
		lambda = struct('t',0,'n',0);

		% TODO make scaling in 1d and 2d interpolation in UnstructuredMesh identical
		% and rename to lambdan and lambdat
		%lambda2 = [0 0]; %1e-4*[1 obj.dwidth];
		%lambda1 = 0;

		%
		% input data
		%

		% error estimation mode (none, standard, jacknife)
		errmode = 'none';

		% externally provided bottom profile
		external_bottom;

		topofbank;
		nprefilter = 0;
	
	%
	%	internal constants
	%
		
		% velocity grid width along cross section (N)
		% TODO, should be dn
		dw;
		% for filtering of shear velocity
		% TODO, this is related to lambda
		% TODO make this part of adcp and filter before binning
		dwf = 0;
		% velocity grid cell height (Z)
		dz;

		% time step (grid cell size in along time vector)
		dt;

		% cross section definition (USER or AUTOMATIC)
		cmethod  = 'USER';

		T_fourier

		n_window = struct('backscatter',1);

		transect;

		%extrapolate = true;
	end % properties
	methods (Static)
		[vval, id]      = extrapolate_velocity_2d_STATIC(vval,bval,cX2,cXX,cYY,mask,vmethod);
		U               = extrapolate_velocity_1d_STATIC(N,U,zb,zs);
		[Bs, id_]       = extrapolate_backscatter_2d_STATIC(Bs,bval,cX2,cXX,cYY,mask,clip);
		[Q, Q_, t0, cdx, msk] = compare(cs, adcp, t0);
		[optstr]        = optstr(obj,varargin);
		%T = split_transect_STATIC(N,w);
		[tid, first, last] = split_transect_STATIC(N,w);
	end % methods Static
	methods
		function obj = CrossSection(varargin)
			obj.grid_n = Grid1();
			if (nargin() < 1)
				% default constructor
				return;
			end

			% take over user defined arguments
			for idx=1:length(varargin)/2
			    field = lower(varargin{2*idx-1});
			    val   = varargin{2*idx};
			    switch (lower(field))
			    case {'xpmax'}
				% this is for backward compatibility
					obj.xlim(2) = val;
			    case {'xpmin'}
					obj.xlim(1) = val;
			    case {'ypmax'}
					obj.ylim(2) = val;
			    case {'ypmin'}
					obj.ylim(1) = val;
			    case {'xlim','ylim'}
					obj.(field) = val;
			    case {'t_max'}
					obj.T_max = val;
			 %   case {'clip'}
			%		% TODO, clip option is deprecated
			%		obj.clip = val;
			    case {'t_rel'}
					obj.T_rel = val;
			    case {'bmethod' ...
				  , 'bmethod2' ...
				  , 'vdimension' ...
				  , 'bfunc' ...
				  , 'vfunc' ...
				  , 'pfunc' ...
				  , 'tmode' ...
				  , 'level' ...
				  , 'dt' ...
				  , 'topofbank'...
				  , 'T_fourier' ...
				  }
				  % ...
				  %'bankful_offset'}
					obj.(field) = val;
				case {'grid_n'}
					obj.grid_n = objcopy(val);
			    otherwise
				obj = setfield_deep(obj,varargin{2*idx-1},val);
			    end % switch field
			end % for idx
		end % constructor
		
		%
		% pseudo member variables, scalars
		%
		% TODO, should be a function of transect
		% TODO also return T, distance and pass x and y
		function [N, obj] = N(obj)
			switch (class(obj.grid_n)) %class(obj.grid_n))
			case {'SparseMesh1'}
				%N = obj.grid_n.val.
				nn = round(obj.transect.dwidth/obj.dw)-1;
				%N  = obj.transect.dwidth*innerspace(-((1:nn)'/(nn+1)-0.5)
				N = obj.transect.dwidth*innerspace(-1/2,1/2,nn);
			case {'Grid1'}
				N = obj.grid_n.cX1;
			case {'RegularizedInterpolator1'}
				N = obj.grid_n.mesh.X;
			otherwise
				error('here');
			end
		end

		% TODO, accept T
		function [X, Y, obj] = XY(obj,N,T)
			if (nargin()<2)
				N = obj.N;
			end
			c   = obj.transect.centre;
			dir = obj.transect.dir; 
			X = c(1) - N*dir(1);
			Y = c(2) - N*dir(2);
		end

		%
		% across section, depth averaged, time independent
		%
		function [zb, obj] = zb(obj,N)
			% TODO pass X and interpolate
			switch (class(obj.grid_n))
			case {'SparseMesh1'}
				zb = obj.var_n('zb');
			case {'Grid1'}
				if (nargin()<2)
					zb = obj.grid_n.val.zb(:,1);
				else
					zb = interp1(obj.grid_n.cX1,obj.grid_n.val.zb(:,1),N,'linear');
				end
			case {'RegularizedInterpolator1'}
				zb = obj.grid_n.vali.zb;
			otherwise
				error('here')
			end
		end % zb

		%
		% cross sectionally averaged, time dependent
		%
		
		% cross sectionally averaged velocity at t
		function [val, obj] = U_t(obj,varargin)
			% TODO also for V,W
			val = cvec(obj.discharge(varargin{:}))./cvec(obj.area(varargin{:}));
		end

		function [val_t, val_tn, val_tnz, obj] = bs_t(obj,varargin)
			% TODO also for V,W
			[val_t] = cvec(obj.var_t('bs_integrated',varargin{:}))./cvec(obj.area(varargin{:}));
			val_tn  = [];
			val_tnz = [];
		end

		function [val_t, val_tn, val_tnz, obj] = Qbs_t(obj,varargin)
			% TODO also for V,W
			[val_t, val_tn] = obj.var_t('qbs',varargin{:}); %./cvec(obj.(varargin{:}));
			val_tnz = [];
		end

		function [val_t, val_tn, val_tnz, obj] = Qss_t(obj,varargin)
			% TODO also for V,W,bs
			[val_t, val_tn] = obj.var_t('sediment_discharge',varargin{:}); %./cvec(obj.(varargin{:}));
			val_tnz = [];
		end


		% cs area at time t
		function [area, obj] = area(obj,t)
			if (nargin()<2)
				t = obj.t0;
			end
			zs = obj.level_t(t);
			switch (class(obj.grid_n))
			case {'SparseMesh1'}
				% TODO water level
				zb   = obj.var_n('zb');
				area = -obj.transect.dwidth*mean(zb);
				area = area*ones(1,length(t));
			case {'0d'}
				valid = obj.transect.valid;                     
		                area = interp1(obj.transect.time(valid),obj.area.total(valid),t,'linear');
			case {'Grid1', 'RegularizedInterpolator1'}
				zb = obj.zb;
				% TODO dangerous
				dn = obj.grid_n.dx1;
				area = csarea(zs,zb,dn);
			otherwise
				error('unimplemented vmethod');
			end % switch
		end % area

		function [perimeter, obj] = perimeter(obj,t)
			if (nargin()<2)
				zs = obj.level_t(obj.t0);
			else
				zs = obj.level_t(t);
			end
			zb = obj.zb;
			% TODO dangerous
			dn = obj.grid_n.dx1;
			perimeter = csperimeter(zs,zb,dn);
		end % perimeter

		function [width, obj] = width(obj,t)
			if (nargin()<2)
				zs = obj.level_t(obj.t0);
			else
				zs = obj.level_t(t);
			end
			zb     = obj.zb;
			dn     = obj.grid_n.dx1;
			width = cswidth(zs,zb,dn);
		end % width

		function [radius, obj] = radius(obj,t)
			if (nargin()<2)
				zs = obj.level_t(obj.t0);
			else
				zs = obj.level_t(t);
			end
			zb     = obj.zb;
			dn     = obj.grid_n.dx1;
			radius = csradius(zs,zb,dn);
		end % radius


		function [lt, obj] = level_t(obj, t)
			if (isempty(obj.level) || isempty(obj.level.time))
				lt = zeros(size(t));
			else
			fdx = isfinite(cvec(obj.level.time)) & isfinite(cvec(obj.level.val));
			switch (sum(fdx))
			case {0}
				error('');
			case {1}
				lt = obj.velocity.val(fdx);
			otherwise
				if (nargin()<2)
					t = obj.t0;
				end
				%lt  = interp1(obj.level.time(fdx), obj.level.val(fdx),t,'linear',NaN);
				dt_max = 0.1*range(obj.t0);
				lt     = interp1_limited(obj.level.time(fdx), obj.level.val(fdx), t, dt_max, true, 'linear','extrap');
			end % switch
			end % if
		end % level_t


		%function [t, n, z, obj] = sigma2z(obj,t,n,s)
		function [z, obj] = sigma2z(obj,t,n,s)
			% TODO scale and add water level
			[A, fdx]  = obj.grid_n.mesh.interpolation_matrix_1d(n);
			zb      = NaN(size(s));
			zb(fdx) = A*obj.grid_n.vali.zb;
			z       = zb.*(1.0-s);  
			z(~fdx) = 0;
		end

		%
		% short-hand functions for tn (depth averaged)
		% varargin can be t

		% spanwise velocity
		function [val, obj] = U_tn(obj,varargin)
			val = obj.var_tn('U',varargin{:});
		end
		% spanwise velocity
		function [val, obj] = V_tn(obj,varargin)
			val = obj.var_tn('V',varargin{:});
		end
		% vertical velocity
		function [val, obj] = W_tn(obj,varargin)
			val = obj.var_tn('W',varargin{:});
		end
		% specific discharge
		function [val, obj] = q_tn(obj,varargin)
			val = obj.var_tn('Q',varargin{:});
		end
		% rouse parameter
		function [val, obj] = rouse_tn(obj,varargin)
			val = obj.var_tn('rouse',varargin{:});
		end
		% backscatter
		function [val, obj] = bs_tn(obj,varargin)
			%val = obj.var_tn('backscatter',varargin{:});
			val = obj.var_tn('bs_integrated',varargin{:});
		end
		% backscatter flux
		function [val, obj] = qbs_tn(obj,varargin)
			val = obj.var_tn('qbs',varargin{:});
		end
		% shear velocity
		function [val, obj] = us_tn(obj,varargin)
			val = obj.var_tn('u_s',varargin{:});
		end
		% logarithmic roughness length
		function [val, obj] = ln_z0_tn(obj,varargin)
			val = obj.var_tn('ln_z0',varargin{:});
		end
		% wake parameter
		function [val, obj] = wake_tn(obj,varargin)
			%val = obj.var_tn('wake',varargin{:});
			val = obj.var_tn('perturbation',varargin{:});
		end

		% local depth
		function [Htn, obj] = d_tn(obj,varargin) 
			l = obj.level_t(varargin{:});
			
			% level  = obj.Ht(t);
			zb     = obj.zb;
			Htn    = bsxfun(@minus,rvec(l),cvec(zb));
		end

		%
		% tnz (time, across, up)
		% 
		
		% short-hand

		% all velocity
		function val_tnz = UVW_tnz(obj,varargin)
			val_tnz = obj.var_tnz('UVW',varargin{:});
		end
		% streamwise velocity
		function val_tnz = U_tnz(obj,varargin)
			val_tnz = obj.var_tnz('Utnz',varargin{:});
		end
		% spanwise velocity
		function val_tnz = V_tnz(obj,varargin)
			val_tnz = obj.var_tnz('Vtnz',varargin{:});
		end
		% vertical velocity
		function val_tnz = W_tnz(obj,varargin)
			val_tnz = obj.var_tnz('Wtnz',varargin{:});
		end
		% backscatter
		function val_tnz = bs_tnz(obj,varargin)
			val_tnz = obj.var_tnz('backscatter_tnz',varargin{:});
		end

		% depth of cell centres
		function [DD, obj] = Dtnz(obj,t)
			switch (class(obj.grid_n))
			case {'0d'}
				error('unavailable');
			case {'Grid1'}
				switch (lower(obj.vdimension))
				case {'n','tn'}
					error('unavailable');
				case {'nz','tnz'}
					% level at time t
					dh    = obj.Ht(t);
					% mean depth of mesh cell centres
					DD0   = obj.grid_nz.cXX2;
					% mean depth
					H0    = obj.zb;
					% scale from mean depth to instantaneous depth
					scale = (H0+dh)./H0;
					% instantaneous depth of mesh cell centres
					DD    = bsxfun(@times,DD0,scale);
				otherwise
					error('here');
				end
			case {'RegularizedInterpolator1'}
				error('not yet implemented');
			otherwise
				error('invalid vmethod');
			end
		end % Dtnz

		%
		% T (depth and cross section averaged)
		%
		% mean depth
%		function Ht = Ht(obj,t)
%			Ht  = interp1(obj.adcp.time,obj.ens.level,t);
%			Ht  = interp1(obj.level.time,obj.level.val,t,'linear');
%		end
%		function [U obj] = U(obj)
%			if (strcmp('NOTIME',upper(obj.cs(1).tmode)))
%				U        = obj.cs.discharge.total(1)./obj.cs.area.total;
%			else
%				U        = obj.cs.discharge.total./obj.cs.area.total;
%			end
%		end
		%function [V obj] = V(obj)
		%	if (strcmp('NOTIME',upper(obj.tmode)))
		%		V        = cvec(obj.discharge.total(2))./obj.area.total;
		%	else
		%		V        = cvec(obj.discharge.V.total)./obj.area.total;
		%	end
		%	% redefinition of the discharge
		%	%obj.cs.discharge.total(2)/obj.cs.area.total(1);
		%end
		%function [W obj] = W(obj)
		%	W = [];
		%	%W        = obj.cs.discharge.total(3)/obj.cs.area.total(1);
		%end


		% turbulent (effective) reynolds number
		%function [Re_t obj] = Re_t(obj)
		%	Re_t = 13*obj.Chezy(1)/Constant.g;
		%end
	end % methods
end % class CrossSection
		
