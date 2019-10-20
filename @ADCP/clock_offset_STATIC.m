% Do 14. Mai 15:26:38 CEST 2015
% Sun 30 Jul 11:51:37 CEST 2017
%% dt : median difference between adcp clock and UTC
%% sd_dt : standard error of dt
% TODO merge with verify_pc_time
function [dt, sd_dt, time, t_gps] = clock_offset_STATIC(adcp,dt0)
	if (nargin() < 2)
		dt0 = [];
	end
	%id = unique(adcp.FileNumber);
	m = max(adcp.FileNumber);
	dt = NaN(m,1);
	sd_dt = NaN(m,1);
	% time of instrument clock (if not overwritten by WinRiver)
	time  = datenum(adcp.timeV);
	% GPS time (UTC)
	date  = adcp.nFiles.ZDA.date;
	t_gps = adcp.nFiles.GGA.UTCtime*[1/24, 1/1440, 1/86400]' + datenum(date);

%	v=vadcp.nFiles.ZDA.UTCtime;
%	d=vadcp.nFiles.ZDA.date;
%	t_gps = datenum(d(:,1),d(:,2),d(:,3),v(:,1),v(:,2),v(:,3));

	for idx=1:m
		fdx = adcp.FileNumber == idx;
		% apply clock offset correction
		if (~isempty(dt0))
			time(fdx) = time(fdx)-dt0(idx);
		end
		if (any(fdx))
			[dt(idx), sd_dt(idx)] = median_man(time(fdx)-t_gps(fdx));
		end
	end
end

