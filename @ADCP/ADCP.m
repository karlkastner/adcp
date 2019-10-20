%% ADCP superclass
%% converts ADCP fixed integer raw data to floats with SI units
%% provides functions for ADCP data manipulation
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
classdef ADCP < handle
	properties (Constant)
		% slope of depth range (m)
		DEPTH_SLOPE     = 0.01;

		% error value of the raw velocity data (SHRT_MIN)
		ERRVAL_VEL    = int16(-32768); 

		% raw attitude angle slopes
		HEADING_SLOPE = 1e-2;
		% to distinguihs instrument from cell depth, instrument depth starts with i
		% TODO check this
		IDEPTH_SLOPE_M_PER_BAR = 100/Constant.g;
		% offset for instrument depth
		% the instrument has already be converted for atm pressure
		IDEPTH_OFFSET_M = 0.0;
		PITCH_SLOPE   = 1e-2;
		% raw pressure slope
		PRESSURE_SLOPE_BAR = 1e-4;
		% raw pressiure offset (this is zero, as the instrument is already corrected for atm pressure)
		PRESSURE_OFFSET_BAR = 0.0;
		% threshold for erronous pressure values (this basically happens when the sensor is in air)
		PRESSURE_ERR_THRESH = 1e9;
		ROLL_SLOPE    = 1e-2;

		% slope of temperature (°K)
		TEMPSLOPE_K  = 0.01;

		% Temperature coefficients
		TEMP_ADC = [9.82697464e1, -5.86074151382e-3, 1.60433886495e-7, -2.32924716883e-12];

		% slope of the raw velocity data (m/s)
		VEL_SLOPE_M_PER_SEC = 1e-3;
		HEADING_MAX = 35999;

		% upper limit for declaring velocity as invalid
		VMAX = 10;
	end % properties (Constant)
	properties (Abstract, Constant)
		% beamangle (1,1) angle between vertical and ADCP beam
		DEFAULT_BEAMANGLE_RAD;
		% slope of echo intensity	
		DEFAULT_KC_DB;
		% first offset of echo intensity
		DEFAULT_C_DB;
		% noise level (second offset of echo intensity)
		% TODO, this is taken from the VADCP is this correct ?
		DEFAULT_EN_DB;
		% operation frequency (Hz)
		FREQ_HZ;
		% tranceiver radius (m)
		% note: the red membrane has a radius of 7cm, but the transducer radisu is only 5.5cm
		AT_M;
		% number of beams
		nbeams;

	end
	% properties shared by all ADCPs
	properties
		% dat : (struct) ADCP raw data as obtained from the adcptools script
		dat;

		%
		% scalar quantities
		%

		% Kc  : (1,1)   slope of echo intensity E in dB (0.45 dB)
		Kc;

		% C   : (1,1)   intercept (offset) of E in dB (-129.1 dB)
	    	C;

		% En  : (1,1) noise level  (offset of E in unit 1)
		En;

		% either HADCP, VADCP or SPADCP
		type;
		
		%
		% bin quantities  [nbin x nens x nbeam]
		%
		% velocity
		velocity = struct('beam',[],'inst',[],'ship',[],'earth',[]);

		% object of bin class
		bin;

		% echo intensity
		echo;

		% TODO, duplicate from utm
		% utm X coordinate
		X;

		% utm Y coordinate
		Y;

		% transversal coordinate of ensemble
		N;

		% streamwise coordinate of ensemble
		T;

		% position of beam intersection with the bottom
		% TODO, should go to ens
		N4;
		T4;

		% depth as measured by each beam
		%H4;

		%
		% beam quantities [nens x nbeam]
		%
		% bottom track velocity
		btvel
		gpsvel;

		% distance to bottom
		btrange_
	
		%
		% ensemble quantities
		%
		% TODO should be bin.n and ens.n
		ens

		% vertical distance between pressure sensor and transducer
		id0 = 0;

		% time : (nt, 1) time of ensemble masurement
		% (days since year 0)
		time;

		% wgs84 position
		wgs84;
		utm = struct( 'X', [], 'Y', [], 'zone', []);

		% transducer depth
		d_transducer = 0;

		% index of ensembles with valid velocity data
		% valid velocity values (no error value / out of range)
		% TODO masks should be properties of bin
		velmask;

		% range mask
		rangemask_;

		% firmware fix information
		fix = struct();

		% dir
		dir

		% filter window length (s)
		filter = struct('t', 600 ...
			... % flag to select the echo intensity filtration mode
			... % log, lin, []
			, 'mode', [] ... 
			, 'win', [] ... 
			, 'n', []);

		% model (1,1) number used for backscatter inversion
		%model;
			
		% flag to witch near field correction on or off
		psiflag = 1;

		% selected noise compensation
		noiseflag = 'log';

		% backscatter_log10 : (nr,nt) : ADCP backscatter in geometric space without correcting for sediment attenuation
		backscatter_log10;

		% beta : (nr,nt) : backscatter in arithmetic space
		backscatter;

		dt_max_gap = 10/86400;

		% backscatter object stores calibration data for conversion to sediment concentration
		backscatter_obj
		
		is_facing_upward = false;
	end	% properties

	methods (Static)
		[dt, t, time, t_gps]	   = clock_offset_STATIC(adcp,dt);
	        adcp                       = squeeze_STATIC(adcp,dt,mode,mask);
	        adcp                       = sort_STATIC(adcp);
		serial                     = convert_raw_serial_STATIC(serial);
		time                       = convert_raw_time_STATIC(timeV);
		[vel, btvel]               = convert_raw_velocity_STATIC(vel,btvel);
		[vel, btvel]               = instrument_to_ship_STATIC(vel,btvel,pitch,roll, iflag);
		[Kc, En]                   = load_RSSI_values_STATIC(serial);
		[vel, btvel]               = ship_to_earth_STATIC(vel, btvel, heading, iflag);
		[VEL] = rotate_velocity(d,VEL,iflag);
	end % methods (STATIC)

	methods
		function obj = ADCP(dat,varargin)
 			obj.backscatter_obj = Backscatter();
			if (nargin() > 0)
			% non-empty constructor

			% set optional parameters
			% TODO, privatise
			for idx=1:length(varargin)/2
				obj = setfield_deep(obj,varargin{1+2*(idx-1)},varargin{2*idx});
			end
	
			% take over raw data from adcp-toolbox
			obj.dat = dat;
	
			%
			% conversion of raw data to physical (SI) units
			%
			% obj.assign_file();
			obj.convert_raw_nmea();
			% obj.convert_raw_serial();
			obj.convert_raw_time();
			obj.convert_raw_velocity();
	
			obj.echo = single(obj.dat.ECHO);
			end % if
	
			% ensemble and bin properties
			% must come after raw conversion
			obj.ens = Ensemble(obj);
			obj.bin = ADCP_Bin(obj);
		end % constructor

		% angle of side beams with respect to centre beam [rad]
		function beamangle_rad = beamangle_rad(obj)
			beamangle_rad = deg2rad(single(obj.dat.HADCPbeamangle));
		end
	
		% TODO, for HADCPs, this is different for each beam
		% dr : (1,1) : along beam length of a single cell (m)
		function dr = dr(obj)
			dr = obj.binsize./cos(obj.beamangle_rad);
		end
	
		% transmit current
		function current_A  = current_A(obj)
		    	current_A  = obj.adc_current_slope()*single(obj.dat.ADC(:,1));
		end
	
		% transmit voltage
		function voltage_V = voltage_V(obj)
			v_raw = single(obj.dat.ADC(:,2));
			v_raw(0 == v_raw) = NaN;
		    	voltage_V  = obj.adc_voltage_slope()*v_raw;
		end
	
		% transmitt power
		function power_W = power_W(obj)
			power_W = obj.voltage_V.*obj.current_A;
		end
	
		% sound velocity
		% cw : (nt,1) velocity of sound (m/s) in fresh water at negligible depth,
		% calculated from the water temperature end sound frequency

		function cw = sound_velocity(obj)
			t_C = obj.transducer_temperature_C();
			cw = Constant.sound_velocity_water(t_C);
			% - obj.TEMPOFFSET_K);
		end
	
		% sound attenuation
		% attenuation by water (dB/m) in fresh water at negligible depth
		%% calculated from the water temperature and sound frequency
		function alpha_w = sound_absorption_water(obj)
			t_C = obj.transducer_temperature_C();
			alpha_w = Constant.sound_absorption_water(obj.FREQ_HZ, 0, 0, t_C);
		end
	
		% beam spreading
		%	at      = 5.5e-2;
		% centerline to half signal strength (not power?)
		% according to ASTM E-1065
		function beamspreading_rad = beamspreading_rad(obj)
			beamspreading_rad = asin(1.2*obj.cw / (2*obj.AT_M*obj.FREQ_HZ));
		end

		%
		% wrapper for static methods
		%
		function [mask, obj] = mask(obj)
			mask = obj.rangemask & obj.velmask;
		end

		function obj = convert_raw_nmea(obj)
			if (isfield(obj.dat,'NMEAGGA'))
				% TODO no magic numbers (zone)
				[obj.wgs84, obj.utm] = nmea2utm(obj.dat.NMEAGGA,obj.utm.zone);
			elseif (isfield(obj.dat,'GGA'))
				[obj.wgs84, obj.utm] = nmea2utm(obj.dat.GGA,obj.utm.zone);
			elseif (isfield(obj.dat,'dFiles'))
				[obj.wgs84, obj.utm] = nmea2utm(obj.dat.dFiles.GGA,obj.utm.zone);
			else
				warning('No coordinates provided, using bottom track positiong as pseudo coordinates');
				%TODO set btflag
			end
		end % convert_raw_NMEA

		% serial number of the ADCP
		function serial = serial(obj)
			serial = obj.convert_raw_serial_STATIC(obj.dat.serial);
		end

		function obj = convert_raw_time(obj)
			obj.time        = obj.convert_raw_time_STATIC(obj.dat.timeV);
		end

		function obj = instrument_to_ship(obj,varargin)
			[obj.velocity, obj.btvel]	= obj.instrument_to_ship_STATIC(obj.velocity,obj.btvel,obj.pitch_rad,obj.roll_rad, varargin{:});
		end

		function obj = load_RSSI_values(obj)
				[obj.Kc, obj.En] = obj.load_RSSI_values_STATIC(obj.serial);
				if (isnan(obj.Kc))
					obj.Kc = obj.DEFAULT_KC_DB*[1 1 1 1];
				end
				if (isnan(obj.En))
					obj.En = obj.DEFAULT_EN_DB*[1 1 1 1];
				end
		end
		function obj = ship_to_earth(obj,varargin)
			[obj.velocity, obj.btvel]    = obj.ship_to_earth_STATIC(obj.velocity, obj.btvel, obj.heading_rad, varargin{:});
		end
		function [t0, obj] = t0(obj)
			t0 = mean(obj.time);
		end

		% TODO rangemask should be a member of ens
		% last a member of ens
		function [ldx, obj] = last_bin(obj,varargin)
			ldx = sum(obj.rangemask(varargin{:}))';
		end
		function [lastFile, obj] = lastFile(obj)
			lastFile = max(obj.dat.FileNumber); 
		end
		function [nt, obj] = nt(obj)
			nt = length(obj.time);
		end
	end % methods
end % classdef ADCP

