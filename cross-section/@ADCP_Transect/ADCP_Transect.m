% Sun 27 May 11:54:42 CEST 2018
%
%% zero dimensional processing of ADCP data
%% no resampling, meshing or gridding
%
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
classdef ADCP_Transect < handle
	properties
		%vadcp
		tdx

		mode = 'returning';

		% end point coordinates
		xlim
		ylim

		% mid-times of transects
		% t0
		time = [];
		% n    = [];

		% first and last ensemble indices
		first
		last
		first_valid
		last_valid
		isvalid
		rid

		velocity
		discharge
		area
		width
		width_limited
		sig

		% maximum streamwise distance for assignment of ensembles
		T_max = 100;
	end % properties
	methods
		function obj = ADCP_Transect(varargin)
			if (nargin() < 1)
				% default constructor
				return;
			end

			% take over user defined arguments
			for idx=1:length(varargin)/2
			    field = lower(varargin{2*idx-1});
			    val   = varargin{2*idx};
			    switch (lower(field))
			    case {'xlim','ylim'}
					obj.(field) = val;
			    case {'t_max'}
					obj.T_max = val;
%			    case {'t_rel'}
%					obj.T_rel = val;
			    otherwise
				obj = setfield_deep(obj,varargin{2*idx-1},val);
			    end % switch field
			end % for idx
		end % constructor

		function obj=import_mmt(obj,mmtname)
			[obj.first,obj.last] = RDI_mmt.read(mmtname);
			% TODO calc auxilary quantities
		end

		function dir = veldirection(obj)
			% compass angle is to right with respect to north,
			% so indeed atan(x,y)
			dir = atan2(obj.velocity(:,1), obj.velocity(:,2));
		end

		% horizontal velocity component only
		function velmag = velmag(obj)
			% compass angle is to right with respect to north,
			% so indeed atan(x,y)
			velmag = hypot(obj.velocity(:,1), obj.velocity(:,2));
		end

		% number of transects
		function n = n(obj)
			n = length(obj.first);
		end

		function [c, obj] = centre(obj)
			c = [mean(obj.xlim); mean(obj.ylim)];
		end % centre

		% width of the domain
		function [dwidth, obj] = dwidth(obj)
			[dir, dwidth] = obj.dir();
		end % dwidth

		function [dir, dwidth, obj] = dir(obj)
			dir = [diff(vertcat(obj.xlim)');
			       diff(vertcat(obj.ylim)')];
			% normalize
			dwidth = hypot(dir(1,:),dir(2,:));
			dir    = dir./dwidth;
		end % dir

		function nlim = nlim(obj)
			nlim = 0.5*obj.dwidth*[-1,1];
		end
	end % methods
end % classdef

