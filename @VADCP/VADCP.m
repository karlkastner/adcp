% Karl Kästner, Berlin

%% coverts raw data of vertical ADCPs into physical quantities
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
classdef VADCP < ADCP
	properties (Constant)
	% definition of abstract ADCP properties
		% beamangle (1,1) angle between vertical and ADCP beam
		DEFAULT_BEAMANGLE_RAD  = deg2rad(20);
		% slope of echo intensity	
		DEFAULT_KC_DB          = 0.38;
		% first offset of echo intensity
		DEFAULT_C_DB	       = -129.1;
		% noise level (second offset of echo intensity)
		DEFAULT_EN_DB	       = 40;
		% slope of voltage (V)
		% undoccumented by RDI, taken from Shields 2010
		ADC_VOLSLOPE           = 0.253765;
		% slope of current (A)
		ADC_CURSLOPE           = 0.011451;
		% operation frequency (Hz)
		FREQ_HZ	               = 1.2e6;
		% tranceiver radius (m)
		% note: the red membrane has a 3.5cm radius,
		% but the transducer radius is only 2.7cm
		% reference : RDI primer
		AT_M  	               = 2.7e-2;
		% number of beams
		nbeams                 = 4;
	end % properties (Constant)
	properties
		% bottom track
		bt;
	%
	% model and parameters selected by the user
	%
		% d_k (1,1) depth of virtual reference
		d_k;

		% r_k;
		% index of virtual reference depth
		%nk;
	%
	% measurement quantities
	%

		% E  : (nr,nt) echo intensity
		%E;
		% r0 : (1,1) : along beam distance towards first cell (m)
		%r0;
		% d0 : (1,1) : depth projected distance to first cell (m)
		d0;
		%dr;
		% dd : (1,1) : depth projected length of a single cell (m)
		dd;
		% L   : (1,1) : length of the transmit pulse (m)
		%L;

		% betdepth : (nt, 4) depth, projected distance to the bottom for each beam
		btdepth;
		% btn (nt, 1) index of the last valid  depth cell (or first invalid?)
		btn;


		% pitch of instrument
		%pitch;
		% role attitude angle of instrument
		%roll;
	%
	% calibration values
	%
		% M_ref : (nc,1) reference sediment concentration at calibration points
		%M_ref;
		% d_ref : (nc,1) depth of reference concetrations
		%d_ref;
		% t_ref : (nc,1) time of reference concentration
		%t_ref;
		% r_ref : (nc,1) transducer distance of reference concentrations
		%r_ref;
		% space index of reference concentrations
		%nr_ref;
		% time index of reference concentrations
		%nt_ref;
	%
	% quantities calculated from the measured quantities
	%

		% nerr  : (1,1) : l2-norm of the regression error (mg/l)
		nerr;
		% M_reg : (nc,1) : sediment concentration at calibration points
		M_reg;
		% corrMM (1,1) : correlation of reference concentration and concentration obtained by inversion
		corrMM;
		% corred (1,1) : correlation of the inversion error with depth
		corred;
		% corret (1,1) : correlation of the inversion error with time
		corret;
		% err   : (nc,1) error of the regressed values copared to calibration values (mg/l)
		err;
		% rmse  : (1,1) root mean square error of regression (mg/l)
		rmse;
		% chi2r : (1,1) : chi^2 statistic of the error
		chi2r;
		% flag to indicate non-normality of the error
		abnormal;
		% mean, standard deviation and bias of the estimated sediment concentration (mg/l) 
		M = struct('mean',[],'std',[],'bias',[]);
		% param : calibrated parameters obained by regression
		%         physical types of parameters depend on the regression model
		param = struct('mean',[],'cov',[],'corr',[],'std',[]);
		% V : (nr, nt) velocity in cartesian coordiantes referenced to the earth
		V;
		% V : (nr, nt) as V, but velocities but cell below bottom are NaN	
%		Vvalid;
		% vbar : (nt,1) depth averaged velocity
		vbar;
		% vbar_f (nt,1) depth averaged velocity, time filtered
		vbar_f;
		% qw   : (nt,1) discharge
		qw;
		% qw_f : (nt,1) time filtered discharge
		qw_f;
		% qc   : (nt,1) sediment transport
		qc;
		% qc_f : (nt,1) time filtered sediment transport
		qc_f;
		% cbar : (nt,1) depth averaged sediment concentration
		cbar;
		
	%
	% internal variables
	%
		% cell of parameters estimated with jacknife
		% 1 : leave one out, 2 : leave two out, ...
		J;

		% filter weight for backscatter filter
		pf_bs = 0.3;

		% nf : filter window length
		nf;

		filt = 0;

		% f : filter window in time
		f;

		% cell array with regression functions of all models
		regC;
		% cell array with evaluation functions of all models
		funcC;
			%
			%X;
			%
			%PC;
			%PC_ref;
			%
			%g_K;
			% japaram: calibrated parameters improved by Jackknife
			%jparam = struct('mean',[],'cov',[],'corr',[],'std',[]);
			% Jparam (nc,nparam) jackknife matrix,
			%        i-th row holds calibration parameters obtained by ignoring calibration point i
			%Jparam;
			% leave-2 out jackknife matrix
			%Jparam2;

			%
		nFiles;
		ibtflag = false;
		ShipReference = 'bt';
		velocity_profile;
		sediment_concentration_profile;

		shear_stress_field = 'sw';
	end % properties

	methods (Static)
		[vel, btvel]       = beam_to_instrument_STATIC(vel,btvel,beamangle_rad,mode,flag);
		 vel               = correct_for_platform_velocity_STATIC(corstr,vel,btvel,FileNumber,iflag);
		[vel, btvel]       = to_beam_STATIC(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad);
		[vel, btvel]       = to_earth_STATIC(corstr,FileNumber,vel,btvel,heading_rad,pitch_rad,roll_rad);
		[btrange, rmask]   = filter_btrange_STATIC(time,vel,btrange,D,roll,pitch,heading,FilenNumber,beamangle_rad);
		[X, Y, S, dX, dY, dS]  = bottom_track_STATIC(time, btvel);
		[str] = optstr(opt);
		q                  = depth_integrate(q,Z,vel,binmask,H,ensmask);
	end % methods (Static)

    methods
	function obj = calc_gpsvel(obj)
		% TODO this is not necessary if there is a VTG field
		if (~isempty(obj.utm.X))
		t = obj.time;
		X = obj.utm.X;
		Y = obj.utm.Y;
		dx = diff(X);
		dy = diff(Y);
		dt = Constant.SECONDS_PER_DAY*diff(t);
		% negate to have ident sign convention as btvel
		obj.gpsvel.earth = -single([dx./dt dy./dt zeros(length(dt),2)]);
		obj.gpsvel.earth(end+1,:) = NaN;
		end
	end
        %
        % wrappers for static methods
        %
	function echo = echo(obj)
		echo = single(obj.dat.ECHO);
	end
        function obj = beam_to_instrument(obj)
            [obj.velocity, obj.btvel] = obj.beam_to_instrument_STATIC( ...
					obj.velocity, obj.btvel,obj.beamangle);
        end % beam_to_instrument
        function obj = correct_for_platform_velocity(obj)
	     % note : if this function is called twice the correction is undone according to ibtflag
	v1 = obj.velocity.earth;
	     switch (obj.ShipReference)
		case {'bt'}
             		obj.velocity = obj.correct_for_platform_velocity_STATIC( ...
				obj.dat.corstr,obj.velocity,obj.btvel,obj.dat.FileNumber,obj.ibtflag);
	     	case {'vtg'}
          		   obj.velocity = obj.correct_for_platform_velocity_STATIC( ...
				obj.dat.corstr,obj.velocity,obj.nFiles.vtgvel,obj.dat.FileNumber,obj.ibtflag);
		otherwise
			error('VADCP correct_for_platform_velocity');
	     end
	v2 = obj.velocity.earth;
	nannorm(v1(:,:,1)-v2(:,:,1))
	     obj.ibtflag = ~obj.ibtflag;
        end % correc_platform_velocity
        function obj = to_earth(obj)
            [obj.velocity, obj.btvel]   = obj.to_earth_STATIC(obj.dat.corstr, ...
		obj.dat.FileNumber,obj.velocity,obj.btvel,obj.heading_rad, ...
		obj.pitch_rad,obj.roll_rad);
        end % to_earth
        function obj = to_beam(obj)
            [obj.velocity, obj.btvel]   = obj.to_beam_STATIC(obj.dat.corstr, ...
		obj.dat.FileNumber,obj.velocity,obj.btvel,obj.heading_rad, ...
		obj.pitch_rad,obj.roll_rad);
        end % to_earth
%        function obj = filter_btrange(obj)
%		if (isfield(obj.dat,'btrange'))	
%	    [obj.btrange, obj.rangemask] = VADCP.filter_btrange_STATIC( ...
%		obj.time, obj.velocity, ...
%		obj.btrange, ...
%		obj.Dt, ...
%		obj.roll_rad, ...
%		obj.pitch_rad, ...
%		obj.heading_rad, ...
%		obj.dat.FileNumber, ...
%		obj.beamangle_rad(1));
%		end
%        end % filter_btrange
        function obj = bottom_track(obj)
	    if (isfield(obj.dat,'btvel'))
	    for idx=1:obj.dat.FileNumber(end)
	   	 fdx = find(idx == obj.dat.FileNumber);
	            [obj.bt.X(fdx,1), obj.bt.Y(fdx,1), obj.bt.S(fdx,1), ...
			 obj.bt.dX(fdx,1), obj.bt.dY(fdx,1), obj.bt.dS(fdx,1)] ...
                    = obj.bottom_track_STATIC(obj.time(fdx), obj.btvel.earth(fdx,:));
	    end
	    end
	end % bottom_track

	    % constructor
	    function obj = VADCP(varargin)
		% call parent constructor
		obj@ADCP(varargin{:});

		% non-empty constructor
		if (~isempty(varargin))


	    	% backscatter of all four beams
		% TODO make a function
	    	%obj.E = single(obj.dat.ECHO);

		% load Kc from library, if known
		obj.load_RSSI_values();
		obj.C = obj.DEFAULT_C_DB;


		obj.convert_nFiles();
%		obj.sediment_concentration_profile = Rouse_Profile();
%		obj.velocity_profile = Log_profile();

		%
		% postprocess data
		%

		% convert VADCP velocity to earth reference
		% must preceed bottom track
		obj.to_earth();

		% velocity from gps track, if not VTG field present
		obj.calc_gpsvel();

		% must preceed filter velocity, as correction may create invalid values
		% must preceed to_earth as correction is only applied to source reference
	        obj.correct_for_platform_velocity();

		% relative coordinates from bottom velocity
		obj.bottom_track();

		% select coordinate source
		if (~isempty(obj.utm.X))
			obj.X = obj.utm.X;
			obj.Y = obj.utm.Y;
		else
			if (~isempty(obj.bt))
				obj.X = obj.bt.X;
				obj.Y = obj.bt.Y;
			end
		end

		% must precede filter_btrange
		% TODO, make this a member of bin
		%obj.calc_bin_range(); 

		% filter depth below instrument
		%obj.filter_btrange();
		
		% must come after filter_btrange
		obj.ens.calc_beamcoords();

		% compsensate for transdurce depth
		% TODO, should go into pressure function
		%obj.idepth_m = obj.idepth_m + obj.id0;

		% filter XY-coordinates
		% has to precede filter vel, as it interpolates coordinates
			% obj.filter_coordinates();


		% filter bin values
		obj.filter_velocity('earth');

		%
		% backscatter calculation
		%
		% TODO, this filters twice if called again later
		obj.calc_backscatter();

		% choices of backscatter inversion models
		% TODO this should go into properties
		%obj.regC =  { @reg1,  @reg2,  @reg3,  @reg4,  @reg5,  @reg6,  @reg7,  @reg8,  @reg9,  @reg10,  @reg11, @reg12, @reg13, @reg5b,  @reg5c};
		%obj.funcC = {@func1, @func2, @func3, @func4, @func5, @func6, @func7, @func8, @func9, @func10, @func11, @func12, @func13, @func5b, @func5c};
		end 
	    end % constructor VADCP()

	    % average backscatter intensity
	    function beta = beta_mu(obj)
		% beam averaging in arithmetic space
			beta   = mean(obj.beta,3);
	    end
	    % standard deviation of backscatter intensity
	    function beta = beta_std(obj)
		% beam averaging in arithmetic space
			beta   = std(obj.beta,3);
	    end
	end % methods
end % class VADCP

