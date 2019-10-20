% Wed Oct 29 11:01:55 CET 2014
% Karl Kastner, Berlin
%
%% container for ADCP ensemble data and properties
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
classdef Ensemble < handle
	properties
		adcp

		% associated transect (all crossigs of one transect)
		tid

		% associated crossings
		cid

		% was vel
		velocity;

		% specific discharge
		% was q
		discharge;

		%
		sediment_discharge
		qbs

%		qbs;
%		bs_integrated;

%		% coordiantes of beam intersection with the bottom
		X4;
		Y4;

		% roughness length fit parameter
		% TODO rename
		param;

		% was us
		shear_velocity;

		% was ln_z0
		ln_roughness_length;

		level;

		transect;

		n_transect;
	end % properties
	methods

	% constructor
	function obj = Ensemble(adcp)
		obj.adcp = adcp;
		%obj.calc_beamcoords();
	end

	function H4 = H4(obj,msk)
		% TODO, this assumes downward looking deployment
		H4 = obj.Dt4 + obj.adcp.d_transducer;
		if (nargin() > 1)
			H4 = H4(msk,:);
		end
	end

	% Thu Jul 17 19:31:07 WIB 2014
	% Karl Kastner, Berlin
	% converts along-beam distance to bottom into depth
	% corrects for pitch, roll, heading
	function Dt4    = Dt4(obj)
		% TODO allow for changing configuration across files
		btrange = obj.adcp.btrange;
		beamangle_rad = obj.adcp.beamangle_rad;
		Dt4     = obj.adcp.range2depth(obj.adcp.btrange,obj.adcp.roll_rad,obj.adcp.pitch_rad,beamangle_rad(1));
	end % filter_btrange()

	% TODO make this a function
	function l4 = level4(obj,fdx)
		if (nargin()<2)
			l4 = repmat(obj.level,1,4);
		else
			l4 = repmat(obj.level(fdx),1,4);
		end
	end

	% bed level
	function zb = zb(obj,msk)
		zb = obj.level(msk)  - obj.H(msk);
	end

	% bed level for each beam
	function zb4 = zb4(obj,msk)
		zb4  = obj.level4(msk) - obj.H4(msk);
	end

%	function Rb = Rb(obj)
%		R = mean(obj.vadcp.btrange(),2);
%	end

	% auto-expanded pseudo variables
	function [n, obj] = n(obj)
		n = length(obj.adcp.time);
	end
	function [time, obj] = time(obj,msk)
		if (nargin()<2)
			time = obj.adcp.time;
		else
			time = obj.adcp.time(msk);
		end
	end
	function [X, obj] = X(obj,msk)
		if (nargin()<2)
			X = obj.adcp.X;
		else
			X = obj.adcp.X(msk);
		end
	end
	function [Y, obj] = Y(obj,msk)
		if (nargin()<2)
			Y = obj.adcp.Y;
		else
			Y = obj.adcp.Y(msk);
		end
	end
	function [H, obj] = H(obj,msk)
		if (nargin()<2)
			% compute water depth
			H  = nanmean(obj.H4,2);
		else
			H  = nanmean(obj.H4(msk),2);
		end
	end


	function [T, obj] = T(obj,msk)
		if (nargin()<2)
			T = obj.adcp.T;
		else
			T = obj.adcp.T(msk);
		end
	end
	function [N, obj] = N(obj,msk)
		if (nargin()<2)
			N = obj.adcp.N;
		else
			N = obj.adcp.N(msk);
		end
	end
	% last valid bin above the bottom
	function [last, obj] = last(obj,varargin)
		last = obj.adcp.last_bin(varargin{:});	
%		% cos(angle) depth = binsize * id + distmidbin1
%		d0    = obj.adcp.distmidbin1;
%		dd    = obj.adcp.binsize;
%		angle = obj.adcp.beamangle_rad;
%		nbins = obj.adcp.nbins;
%		nens  = obj.n;
%		last  = zeros(nens,1);
%		% for each file
%		for idx=1:obj.adcp.file.n
%			fdx        =  obj.adcp.file.id{idx};
%			%TODO this is wrong, use distance to bottom not H!
%			last(fdx)  = floor( (cos(angle(idx))*obj.adcp.H(fdx) - d0(idx))/dd(idx) );
%			% if depth is deeper than last bin, limit
%			last(fdx)  = max(1,min(nbins(idx),last(fdx)));
%		end
	end % last
	end % methods
end % class Ensemble

