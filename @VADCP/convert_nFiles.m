% Mo 24. Aug 14:38:38 CEST 2015
% Karl Kastner, Berlin
%% convert coordinates of NMEA-nFiles
function obj = convert_nFiles(obj)
	if (isfield(obj.dat,'nFiles'))
		if (isfield(obj.dat.nFiles,'GGA'))
		% hack for south - not properly read by ADCPTOOLS
		% TODO no magic numbers (49M,1e7)
		[X, Y] = latlon2utm(-obj.dat.nFiles.GGA.lat,obj.dat.nFiles.GGA.long,'49M');
		Y = 1e7-Y;
		obj.nFiles.X = X;
		obj.nFiles.Y = Y;
		t = obj.time;
		dx = diff(X);
		dy = diff(Y);
		dt = Constant.SECONDS_PER_DAY*diff(t);
		obj.nFiles.gpsvel.earth = -single([dx./dt dy./dt zeros(length(dt),2)]);
		obj.nFiles.gpsvel.earth(end+1,:) = NaN;
		end

		% [X Y] = latlon2utm(vadcp.dat.nFiles.GGA.lat,vadcp.dat.nFiles.GGA.long,'49M'); 
		if (isfield(obj.dat.nFiles,'VTG'))
			heading = deg2rad(obj.dat.nFiles.VTG.TrackDegTrue);
			mag     = obj.dat.nFiles.VTG.SpeedKmH/3.6;
			vx	= sin(heading).*mag;
			vy	= cos(heading).*mag;
			obj.nFiles.vtgvel.earth = -[vx vy zeros(size(vx,1),2)];
		end % if isfield VTG
	end % if isfield nFiles
end % function convert_nFiles

