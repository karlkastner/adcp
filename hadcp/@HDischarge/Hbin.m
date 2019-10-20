% Fri Jan  9 11:24:09 CET 2015
% Karl Kastner, Berlin
%
% bin container class for horizontal ADCP bin values 
%
% TODO compute z0 in 10m intervals
% TODO solve imag ln_z0
% TODO smooth vadcp input
% TODO smooth + extrapolate z_0
% TODO smooth + extrapolate bottom
% TODO tilt correction for vel profile
% todo, reoder arguments nc, nc x nN, nh*nt

% TODO make this a HADCP superclass
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
classdef Hbin < handle
	properties
		nbin
		nens
		% settings
		method
		worder
		% scalar calibration values
		t0
		q0
		nc
		% values per time step
		level1
		offset1
		% calibration values at HADCP bins, averaged over campaigns (necessary?)
		bottom1
		ln_z01
		nh0
		mh0
		% calibration values at HADCP bin, not averaged over campaigns
		ln_z0_A;
		h0_A
		% TODO : rename into cs.qs0_A
		qs0
		% calibration values over cross section
		% TODO, make cs a struct
		csh0_A;
		csln_z0_A;
		% matrix calibration values
		vel		
%		ln_z0_A

		% calculated values
		cw 
		% expanded values
		pitchoffset
		q_tot
		Q
		Q_err		
	end % properties
	methods

	% Mon Jul 28 15:31:54 WIB 2014
	% Karl Kastner, Berlin
	function obj = Hbin( ...
		t0, q0, Nz, csln_z0_A, Nb, bottom_A, Nq, qs_A, htime, Nh, ...
		level, offset, vel, r, pitch_rad, method, worder)
			%nbin,nens)
			%obj.nbin = nbin;
			%obj.nens = nens;
		%end

		% number of calibration campaigns
		obj.nc      = length(t0);
		% number of hadcp bins alogn the beam
		obj.nbin    = size(vel,1);
		% number of hadcp ensembles in time
		obj.nens    = size(vel,2);
		obj.method = method;
		obj.worder = worder;
		obj.level1  = level;
		obj.offset1 = offset;
		obj.q0 = q0;
		obj.csln_z0_A = csln_z0_A;
	
		% average  roughness length over the calibration campaigns
		% TODO warn if there are NaN values
	%	ln_z0  = nanmean(ln_z0_A,2);
		% TODO this is more or less a quick fix
	%	csln_z0_A(:) = nanmean(csln_z0_A(:));
	
		% interpolate roughness length to hadcp bins
		obj.ln_z0_A = interp1(Nz, csln_z0_A, Nh, 'linear');
		obj.ln_z01   = nanmean(obj.ln_z0_A,2);
	
		% water level at calibration campaigns
		level0  = interp1(htime, level, t0,'linear');
		
		% interpolate bottom profile to hadcp bins
		bottom  = interp1(Nb, bottom_A, Nh, 'linear');
	
		% water depth during calibration campaigns at hadcp bins
		for idx=1:obj.nc
			obj.h0_A(:,idx) = bottom(:,idx) + level0(idx);
		end
	
		% average bathymetry at hadcp bins below level
		obj.bottom1  = mean(bottom,2);
	
		% interpolate specific discharge to hadcp bins
		qs0v	= interp1(Nq, qs_A, Nh, 'linear');
		obj.csh0_A  = bottom_A   + repmat(level0',size(bottom_A,1),1);
	
		% velocity at HADCP bin centre
		obj.velocity = vel;
	
		% pitch offset
		obj.pitchoffset = repmat(r(:),1,obj.nens).*repmat(sin(pitch_rad(:)'),obj.nbin,1);

		% normalised depth
		% this does not effect the discharge estimate, but makes the 
		% coefficients of intersect and slope comparable
		obj.mh0 = mean(obj.h0_A,2);
		obj.nh0 = max(abs(obj.h0_A-repmat(obj.mh0,1,obj.nc)),[],2);


		
		% enpand N coordinate for each ensemble
		%bin.N       = repmat(Nh,1,nens);
	
		% relative depth of the bin centres
		%bin.S       = bin.idepth ./ bin.bdepth;
	
		%bin.c_u   = bin.h.*c_; %/(log(bin.z)-bin.ln_z0).*(log(bin.h)-bin.ln_z0-1)
	
		%vel./(log(bin.z)-bin.ln_z0).*(log(bin.h)-bin.ln_z0-1);
	
		% specific discharges during calibration
		qs0h        = interp1(htime,obj.qs',t0,'linear')';

		% qs0v
		obj.qs0 = qs0h;
	
	end % constructor



		% expand reference level to all hadcp bins
		function [level obj] = level(obj)
			level  = repmat(single(obj.level1(:)'),obj.nbin,1);
		end

		% expand offset to all hadcp bins
		% TODO tilt correction
		function [offset obj] = offset(obj)
			offset  = repmat(single(obj.offset1),obj.nbin,1);
		end

		% expand bathymetry to all hadcp ensembles
		function [bottom obj] = bottom(obj)
			bottom  = repmat(obj.bottom1,1,obj.nens);
		end

		% expand roughness length to all each ensembles
		function [ln_z0 obj] = ln_z0(obj)
			ln_z0 = repmat(obj.ln_z01,1,obj.nens);
		end

		% local water depth
		function[h obj] = h(obj)
			h   = obj.bottom + obj.level;
		end

		% specific HADCP discharge
		function [qs obj] = qs(obj)
			qs     = obj.c_u.*obj.velocity;
		end

		% stage dependend discharge weight
		function [c_q obj] = c_q(obj)
			hi = obj.hi;
			c_q = repmat(obj.cw(:,1),1,obj.nens);
			for jdx=1:obj.worder
				c_q = c_q + repmat(obj.cw(:,jdx+1),1,obj.nens).*hi.^jdx;
			end
		end

		function [c_u obj] = c_u(obj)
			c_u   = obj.h.*obj.c_;
		end

		function [Ubar obj] = Ubar(obj)
			% depth averaged velocity from HADCP
			Ubar   = obj.c_.*obj.velocity;
		end

		% elevation of hadcp bin centres above the bottom
		function [z obj] = z(obj)
			z       = obj.bottom - obj.offset - obj.pitchoffset;
			fdx = (z <= 0);
			z(fdx) = NaN;
		end
		% water level scaled for discharge computation
		function [hi obj] = hi(obj)
			hi = (obj.h-repmat(obj.mh0,1,obj.nens))./repmat(obj.nh0,1,obj.nens);
		end

		% scale velocity to specific discharge
		% TODO make this c_u and replace c_u by c_u h every where else
		function [c_ obj] = c_(obj)
			h = obj.h;
			z = obj.z;
			ln_z0 = obj.ln_z0;
			c_ = (log(h)-ln_z0-1)./(log(z)-ln_z0);
		end
	end % methods
end % class Hbin

