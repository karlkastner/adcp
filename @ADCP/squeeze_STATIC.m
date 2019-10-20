% 2014-09-04 15:56:35.261107329 +0200
% Karl Kastner, Berlin
%
%% cut ensembles, skip ensembles or average ensembles in time
%%
%% adcp : adcp structure
%% dt   : time between output ensembles in seconds
%% mode : {'average', 'skip'}
%% mask : selection of ensembles to keep (computed from dt if not provided)
% TODO handle NMEAGGA for mean
% TODO treat nFiles structure
function out = squeeze_STATIC(adcp,dt,mode,mask)
	dtprint = 10;
	% convert dt to days
	% furthermore the threshold is set to 0.9 times the interval
	% so not to average samples that are almost ident
	dt = 0.9*dt / Constant.SECONDS_PER_DAY;
	% allocate memory
	out   = struct();
	itime = datenum(adcp.timeV);
	nt    = length(itime);
	otime = zeros(nt,1);
	field1 = {};
	field2 = {};
	field3 = {};
	nmea = false;
	nFiles = false;
	for field=rvec(fieldnames(adcp))
		field = field{1};
		if (isstruct(adcp.(field)))
		% NMEA data is a bit more complicated
		if (strcmp(field,'NMEAGGA'))
			out.NMEAGGA.Lat  = zeros(nt,1);
			out.NMEAGGA.Long = zeros(nt,1); 
			% apply negative sign for southern hemispehere, thus all squeezed
			% samples will be referenced with respect to north
			out.NMEAGGA.SN   = repmat('N',nt,1);
			hemisphere = -(adcp.NMEAGGA.SN == 'S') + (adcp.NMEAGGA.SN ~= 'S');
			lat        = nanmean(hemisphere.*adcp.NMEAGGA.Lat,2);
			lon        = nanmean(adcp.NMEAGGA.Long,2);
			nmea = true;
		elseif(strcmp(field,'nFiles'))
			for f1=rvec(fieldnames(adcp.nFiles))
				for f2=rvec(fieldnames(adcp.nFiles.(f1{1})))
					out.(field).(f1{1}).(f2{1}) = zeros(nt,1);
				end
			end
			nFiles = true;
		else
			warning(['discarding input field ' field ' which is of type struct']);
		end % if strcmp
		elseif (isempty(adcp.(field)) || isscalar(adcp.(field)))
			out.(field) = adcp.(field);
		elseif (isvector(adcp.(field)))
			if (length(adcp.(field)) == nt)
				out.(field) = alloc(size(adcp.(field)), class(adcp.(field)));
				field1{end+1} = field;
			else
				out.(field) = adcp.(field);
			end
		elseif (2 == ndims(adcp.(field)))
			if (size(adcp.(field),1) == nt)
				out.(field) = alloc(size(adcp.(field)), class(adcp.(field)));
				field2{end+1} = field;
			else
				out.(field) = adcp.(field);
			end
		else % assume 3d
			if (size(adcp.(field),2) == nt)
				out.(field) = alloc(size(adcp.(field)), class(adcp.(field)));
				field3{end+1} = field;
			else
				out.(field) = adcp.(field);
			end
		end

	end % for field
	
%	jdx = 0;
	% if no mask was provided compute mask based on time step
	if (nargin() < 4 || isempty(mask))
		mask = false(nt,1);
		% last input time
%		t   = itime(1);
		% last input index
		jdx = 1;
		for idx=2:nt
		if ( (itime(idx)-itime(jdx) > dt) ...
                    || (adcp.FileNumber(idx) ~= adcp.FileNumber(jdx)) )
			mask(idx-1) = true;
			jdx = idx;
		end % if
		end % for idx
		% last sample
		mask(end) = true;
	end % if
	
	% reset input index	
	jdx = 0;
	% output index
	odx = 0;
	timer = tic();
	tlast = toc(timer);
	for idx=1:nt
		if ( mask(idx) )
		    %(itime(idx)-t > dt) ...
                    %|| (adcp.FileNumber(idx) ~= adcp.FileNumber(jdx)) )
			% progress status
			if (toc(timer)-tlast > dtprint)
				tlast = toc(timer);
				fprintf(1,'Progress: %g\n%% %gs\n',idx/nt,tlast);
			end % if dtprint exceeded
			odx            = odx+1;
		% squeeze ensemble data
		switch (mode)
		case {'average'}
			otime(odx)     = mean(itime(jdx+1:idx));
			for kdx=1:length(field1)
				out.(field1{kdx})(odx) = mean(adcp.(field1{kdx})(jdx+1:idx));
			end
			for kdx=1:length(field2)
				% TODO, this messses things up in case of error values (0 and inf)
				% this will also give awkward results for bits, strings and number of pings
				out.(field2{kdx})(odx,:) = mean(adcp.(field2{kdx})(jdx+1:idx,:),1);
			end
			for kdx=1:length(field3)
				out.(field3{kdx})(:,odx,:) = mean(adcp.(field3{kdx})(:,jdx+1:idx,:),2);
			end
		case {'skip'}
			otime(odx)     = itime(idx);
			for kdx=1:length(field1)
				out.(field1{kdx})(odx) = adcp.(field1{kdx})(idx);
			end
			for kdx=1:length(field2)
				out.(field2{kdx})(odx,:) = adcp.(field2{kdx})(idx,:);
			end
			for kdx=1:length(field3)
				out.(field3{kdx})(:,odx,:) = adcp.(field3{kdx})(:,idx,:);
			end
		otherwise
			error('unimplemented squeeze method');
		end % switch mode
		if (nmea)
			out.NMEAGGA.Lat(odx)  = lat(idx);	
			out.NMEAGGA.Long(odx) = lon(idx);	
		end
		if (nFiles)
			for f1=rvec(fields(adcp.nFiles))
				for f2 = rvec(fields(adcp.nFiles.(f1{1})))
					out.nFiles.(f1{1}).(f2{1})(odx) = adcp.nFiles.(f1{1}).(f2{1})(idx);
				end
			end
		end
			jdx = idx;
		end % if dt exceeded
	end % for idx
	% truncate output
	for kdx=1:length(field1)
		out.(field1{kdx}) = out.(field1{kdx})(1:odx);
	end
	for kdx=1:length(field2)
		out.(field2{kdx}) = out.(field2{kdx})(1:odx,:);
	end
	for kdx=1:length(field3)
		out.(field3{kdx}) = out.(field3{kdx})(:,1:odx,:);
	end
	if (nmea)
		out.NMEAGGA.Lat  = out.NMEAGGA.Lat(1:odx);
		out.NMEAGGA.Long = out.NMEAGGA.Long(1:odx);
		out.NMEAGGA.SN   = out.NMEAGGA.SN(1:odx);
	end
	if (nFiles)
		for f1=rvec(fields(adcp.nFiles))
			for f2 =rvec(fields(adcp.nFiles.(f1{1})))
				out.nFiles.(f1{1}).(f2{1}) = out.nFiles.(f1{1}).(f2{1})(1:odx);
			end
		end
	end
	out.timeV      = datevec(otime(1:odx));
end % adcp_squeeze

