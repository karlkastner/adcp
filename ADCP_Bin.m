% Sat Jan 31 13:24:00 CET 2015
% Karl Kastner, Berlin
%% ADCP bin (single velocity values)
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
classdef ADCP_Bin < handle
	properties
		adcp
		sediment_concentration
	end % properties
	methods
	function obj = ADCP_Bin(adcp)
		obj.adcp = adcp;
	end

	%
	% auto-expanding variables to save memory
	%
	function [n, obj]   = n(obj)
		n = max(obj.adcp.nbins);
	end

	function [mask, obj] = mask(obj)
		mask = obj.adcp.mask;
	end

	% TODO, this is redundant to rangemask
	function [zmask, obj] = zmask(obj,skip)
		% TODO, skip should be class member
		if (nargin<2)
			skip = 0;
		end
		zmask = true(obj.n,obj.adcp.ens.n);
		%obj.adcp.ens.n,1);
		last = obj.adcp.ens.last;
		last = max(1,(last+1-skip));
		for idx=1:obj.adcp.ens.n
			zmask(last:end,idx) = false;
		end
	end

	function [ret, obj] = time(obj,varargin)
		ret = repmat(rvec(obj.adcp.time(varargin{:})),obj.n,1);
	end

	function [ret, obj] = X(obj)
		% TODO awkward construct, make this sub-classes of both HADCP and VADCP
		if (~isa(obj.adcp,'HADCP'))
			ret = repmat(rvec(obj.adcp.X),obj.n,1);
		else
			[X,Y] = obj.adcp.calc_bin_coordinates();
			ret = X;
		end
	end
	function [ret, obj] = Y(obj)
		%ret = repmat(rvec(obj.adcp.Y),obj.n,1);
		% TODO awkward construct, make this sub-classes of both HADCP and VADCP
		if (~isa(obj.adcp,'HADCP'))
			ret = repmat(rvec(obj.adcp.Y),obj.n,1);
		else
			[X,Y] = obj.adcp.calc_bin_coordinates();
			ret = Y;
		end
	end
	function [ret, obj] = N(obj,cdx)
		if (~isa(obj.adcp,'HADCP'))
		if (nargin()<2)
			ret = repmat(rvec(obj.adcp.N),obj.n,1);
		else
			ret = repmat(rvec(obj.adcp.N(:,cdx)),obj.n,1);
		end
		else
		        %[hadcp.N, hadcp.T]         = xy2nt(hadcp.bin.X, hadcp.bin.Y, opt.cs.centre(1), opt.cs.centre(2), opt.cs.dir);
			[N, T] = xy2nt(obj.X, obj.Y, opt.cs.centre(1), opt.cs.centre(2), opt.cs.dir);
			ret = N;
		end
	end
	function [ret, obj] = T(obj,cdx)
		if (~isa(obj.adcp,'HADCP'))
		if (nargin()<2)
			ret = repmat(rvec(obj.adcp.T),obj.n,1);
		else
			ret = repmat(rvec(obj.adcp.T(:,cdx)),obj.n,1);
		end
		else
		        %[hadcp.N, hadcp.T]         = xy2nt(hadcp.bin.X, hadcp.bin.Y, opt.cs.centre(1), opt.cs.centre(2), opt.cs.dir);
			[N, T] = xy2nt(obj.X, obj.Y, opt.cs.centre(1), opt.cs.centre(2), opt.cs.dir);
			ret = N;
		end
	end
	% depth of bin cell below transducer, without tilt correction
	function [Dt, obj] = Dt(obj,msk)
		Dt    = obj.adcp.Dt;
	
		if (nargin()<2)
			fid = obj.adcp.dat.FileNumber;
		else
			fid = obj.adcp.dat.FileNumber(msk);
		end
		Dt = Dt(:,fid);
	end % Dt

	% projected distance between sufrace and cell centre for each beam
	function [Dt4, obj] = Dt4(obj)
		Dt            = obj.adcp.Dt;
		% TODO file individual
		Dt            = Dt(:,1);
		beamangle_rad = obj.adcp.beamangle_rad;
		Dt4 = range2depth(Dt,obj.adcp.roll_rad,obj.adcp.pitch_rad,beamangle_rad(1));
	end

	% project distance between surface and of cell centre
	function [Ds, obj] = Ds(obj,msk)
		Ds = obj.adcp.Ds;

		if (nargin()<2)
			fid = obj.adcp.dat.FileNumber;
		else
			fid = obj.adcp.dat.FileNumber(msk);
		end

		Ds = Ds(:,fid);
	end % Ds

	% range (euclidean (unprojected) distance) from transducer
	function [R, obj] = R(obj,msk)
		R     = obj.adcp.R;

		if (nargin()<2)
			fid = obj.adcp.dat.FileNumber;
		else
			fid = obj.adcp.dat.FileNumber(msk);
		end
	
		R = R(:,fid);
	end % R

	% depth (surface elevation above bottom, includes transducer depth)
	function [ret, obj] = H(obj,msk)
		if (nargin() < 2)
			ret = repmat(rvec(obj.adcp.ens.H),obj.n,1);
		else
			ret = repmat(rvec(obj.adcp.ens.H(msk)),obj.n,1);
		end
	end % H

	% distance of bin centre above bottom
	% TODO rename into Db
	function [Z, obj] = Z(obj,varargin)
		if (~obj.adcp.is_facing_upward)
			Z = obj.H(varargin{:}) - obj.Ds(varargin{:});
		else
			Z = obj.Dt(varargin{:}) + obj.adcp.d_transducer;
		end
	end % Z

	function [velocity, obj] = velocity(obj)
		velocity = obj.adcp.velocity;
	end

	% S-coordinate (normalised depth)
	function [ret, obj] = S(obj,varargin)
		ret = obj.Z./obj.H;
	end

	% vertical distance above bottom for individual bins
	% TODO rename into Db4
	function Z4     = Z4(obj)
		Dt4     = obj.Dt4;
		Dt4_max = obj.adcp.ens.Dt4;
		Z4      = zeros(size(Dt4));
		for idx=1:4
			Z4(:,:,idx) = bsxfun(@minus,Dt4_max(:,idx).',Dt4(:,:,idx));
		end
	end

	% s-coordinate of bins, with respect to beam individual distance to the
	% bottom
	function S4 = S4(obj)
		Z4 = obj.Z4;
		H4 = obj.adcp.ens.H4;
		S4  = zeros(size(Z4));
		for idx=1:4
			S4(:,:,idx) = bsxfun(@times,Z4(:,:,idx),1./H4(:,idx).');
		end
	end

	function [bs, obj] = backscatter(obj)
		bs = obj.adcp.backscatter;
	end

	end % methods
end % classdef Bin

