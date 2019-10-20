% Sat Nov  2 11:34:55 UTC 2013
% Karl Kästner, Berlin
%% coverts raw data of horizontal ADCPs into physical quantities
%% and provides functions for data processing
% TODO : distance R of side beams is incorrect, bins of side beams are larger !!!
% TODO : the distance to the first cell seems not correctly calculated, see RDI mail
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
classdef HADCP < ADCP
	properties (Constant)
		% beamangle (1,1) angle between vertical and ADCP beam
		DEFAULT_BEAMANGLE_RAD = deg2rad(25);
		% slope of echo intensity	
		DEFAULT_KC_DB         = 0.45;
		% first offset of echo intensity
		DEFAULT_C_DB	      = -139.5;
		% noise level (second offset of echo intensity)
		% TODO, this is copied from the VADCP is this correct ?
		DEFAULT_EN_DB	      = 40;
		% operation frequency (Hz)
		FREQ_HZ	              = 614.4e3;
		% tranceiver radius (m)
		% note: the red membrane has a radius of 7cm,
		% but the transducer radius is only 5.05cm
		% reference : RDI primer
		AT_M	              = 5.05e-2;
		% number of beams
		nbeams                = 3;
		% TODO doubplicte of beam_spreading_rad
		SPREADANGLE_RAD = deg2rad(1.4); %deg2rad(1.3/2);
	end % properties (Constant)
	properties
		%psiflag = 1;
		% flag indicating which beam is used for extended backscatter calculations
		%beamflag = 3;
		% TODO ???
		boot;
		% TODO should be named ens.nping
		pingperens;
	end
	properties % variable properties
		vel_unfixed;

	end % properties
	methods (Static)
		[vel, btvel] = beam_to_instrument_STATIC(vel,btvel,beamangle_rad,mode,iflag);
	        [vel, btvel] = to_beam_STATIC(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad);
		[vel, btvel] = to_earth_STATIC(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad,beamangle_rad,mode);
	        [vel, btvel] = to_instrument_STATIC(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad);
	        [vel, btvel] = to_ship_STATIC(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad);
		[vel, fixed] = firmware_fix_STATIC(corstr,FileNumber,firmver,firmrev,vel,btvel,heading_rad,pitch_rad,roll_rad);
		vel = reorder_velocity_STATIC(corstr,FileNumber,vel);
	end % methods (Static)
	methods
		%
		%  wrappers for static methods
		%
		function obj = beam_to_instrument(obj, mode, iflag)
			[obj.velocity obj.btvel] = HADCP.beam_to_instrument_STATIC( ...
				obj.velocity, obj.btvel, obj.beamangle_rad, mode, iflag);
		end % beam_to_instrument
	        function obj = to_beam(obj)
			dat = obj.dat;
			[obj.velocity obj.btvel] = HADCP.to_beam_STATIC( ...
				dat.corstr,dat.FileNumber,obj.velocity,obj.btvel, ...
				obj.heading_rad,obj.pitch_rad,obj.roll_rad);
		end % to_beam
	        function obj = to_instrument(obj)
			dat = obj.dat;
			[obj.velocity obj.btvel] = HADCP.to_instrument_STATIC( ...
				dat.corstr,dat.FileNumber,obj.velocity,obj.btvel, ...
				obj.heading_rad,obj.pitch_rad,obj.roll_rad);
		end % to_beam
	        function obj = to_ship(obj)
			dat = obj.dat;
			[obj.velocity obj.btvel] = HADCP.to_ship_STATIC( ...
				dat.corstr,dat.FileNumber,obj.velocity,obj.btvel, ...
				obj.heading_rad,obj.pitch_rad,obj.roll_rad);
		end % to_beam
	        function obj = to_earth(obj,mode)
			dat = obj.dat;
			if (nargin() < 2)
				mode = [];
			end
			[obj.velocity obj.btvel] = HADCP.to_earth_STATIC( ...
				dat.corstr,dat.FileNumber,obj.velocity,obj.btvel, ...
				obj.heading_rad,obj.pitch_rad,obj.roll_rad, obj.beamangle_rad,mode);
		end % to_beam
		function obj = firmware_fix(obj)
			dat = obj.dat;
			[obj.velocity fixed] = HADCP.firmware_fix_STATIC(...
				dat.corstr, dat.FileNumber, ...
				dat.firmver, dat.firmrev, obj.velocity, obj.btvel, ...
				obj.heading_rad, obj.pitch_rad, obj.roll_rad);
			% quick fix
			for idx=1:length(fixed)
				if (fixed(idx))
					obj.dat.corstr{idx}(1:4) = 'Beam';
				end
			end
		end % firmware_fix
		function obj = reorder_velocity(obj)
			obj.velocity = HADCP.reorder_velocity_STATIC(...
				obj.dat.corstr,obj.dat.FileNumber,obj.velocity);
		end

		% rotation to cs-vel
		function ret = to_cs(obj,dir)
			obj.dir    = dir;
			R          = [dir(2) dir(1); -dir(1) dir(2)];
			vel        = obj.velocity.earth;
			vel_       = reshape(vel(:,:,1:2),[],2);
			vel_       = vel_*R';
			ret        = reshape(vel_,[size(vel,1),size(vel,2),2]);
			ret(:,:,3) = vel(:,:,3);
			obj.velocity.cs = ret;
		end % to_cs

		% constructor
		function obj = HADCP(dat, varargin)
			% call parent constructor
			obj@ADCP(dat,varargin{:});
			% empty constructor
			if (nargin() < 1)
				return;
			end
		

			% load Kc from library, if known
			obj.load_RSSI_values();
			obj.C = obj.DEFAULT_C_DB;

			% apply manual fixes
			fieldname_C = fieldnames(obj.fix);
			for idx=1:length(fieldname_C)
				obj.(fieldname_C{idx}) = obj.fix.(fieldname_C{idx});
			end

%			obj.R    = obj.R*cos(obj.beamangle); <- TODO this is different for the side beams
			obj.load_RSSI_values();

			% expand file-data
			for idx=1:length(obj.dat.pingperens)
				fdx = obj.dat.FileNumber == idx;
				obj.pingperens(fdx) = dat.pingperens(idx);
			end

			%
			% postprocess data
			%

			% correct inconsistent order of error velocity and void fourth beam
			% put error vel into fourth bin
			obj.reorder_velocity();
			% TODO somehow the error velocity is screwed			
			obj.velocity.earth(:,:,4) = 0;
			obj.velocity_unfixed = HADCP.to_beam_STATIC(obj.dat.corstr,obj.dat.FileNumber,obj.velocity,[],obj.heading_rad,obj.pitch_rad,obj.roll_rad);
			% fix the firware bug
			obj.firmware_fix();	
			obj.calc_bin_range();
			%obj.R = obj.D;
			%obj.calc_bin_coordinates();

			% Note : this is only the mask creation, not yet the mask application			
			obj.filter_range();

			% sound velocity
			%obj.c_water   = Constant.sound_velocity_water(obj.temperature_K - obj.TEMPOFFSET_K);

			%			
			% backscatter calculation
			%
			obj.calc_backscatter();
	
			% determine, which coordinate system is used
			% TODO, this is superfluos
			%for idx=1:length(obj.dat.corstr)
			%	obj.isbeam(idx)  = strcmp(obj.dat.corstr{idx}(1:4),'Beam');
			%	obj.isinstr(idx) = strcmp(obj.dat.corstr{idx}(1:4),'Inst');
			%	obj.isearth(idx) = strcmp(obj.dat.corstr{idx}(1:4),'Eart');
			%end

			% calculate velocities in instrument coordinates
	%		obj.calc_vel_instrument();
			%obj.calc_vel_earth();
	%		obj.calc_vel_beam();
		end % constructor HADCP()

	end % methods
end % classdef HADCP

