% Fr 22. Mai 17:14:39 CEST 2015
% Karl Kastner, Berlin
%% sort files by start time
function adcp = sort_STATIC(adcp)

	itime           = cvec(datenum(adcp.timeV));
	nt              = length(itime);

	% start indices of files
	fdx             = [0; cvec(find(diff(double(adcp.FileNumber))~=0))];

	% length of files
	l               = diff([fdx; nt]);

	% sort files by start time
	filetime         = itime(fdx+1);
	[filetime_, sdx] = sort(filetime);

	% new index of individual samples
	nid = NaN(nt,1);
	n   = 0;
	nf  = length(fdx);
	for idx=1:nf
		%fid = itime(fdx(idx));
		% length of current section
		%dn            = fdx(sdx(idx+1))-fdx(sdx(idx));
		%fdx(sdx(idx+1))-1;
		% new indices
		nid(n+1:n+l(sdx(idx))) = fdx(sdx(idx))+1:fdx(sdx(idx))+l(sdx(idx));
		n = n+l(sdx(idx));
	end
%	v = 1:nt;
%	nid = v(nid);
%	nid = nid(nid);
%	t_ = itime(nid);
	nid(nid) = 1:nt;
	
	field1 = {};
	field2 = {};
	field3 = {};
	nmea   = false;
	nFiles = false;
	for field=rvec(fieldnames(adcp))
		field = field{1};
		if (isstruct(adcp.(field)))
		% NMEA and nFiles requires reordering of sub-fields
		if (strcmp(field,'NMEAGGA'))
			adcp.NMEAGGA.Lat(nid,:)  = adcp.NMEAGGA.LAT;
			adcp.NMEAGGA.Long(nid,:) = adcp.NMEAGGA.LON; 
			adcp.NMEAGGA.SN(nid,:)   = adcp.NMEAGGA.SN;
			nmea = true;
		elseif(strcmp(field,'nFiles'))
			for f1=rvec(fieldnames(adcp.nFiles))
				for f2=rvec(fieldnames(adcp.nFiles.(f1{1})))
					adcp.(field).(f1{1}).(f2{1})(nid) = ...
                                                  adcp.(field).(f1{1}).(f2{1});
				end
			end
			nFiles = true;
		else
			warning(['discarding input field ' field ' which is of type struct']);
		end % if strcmp
		elseif (isempty(adcp.(field)) || isscalar(adcp.(field)))
			% nothing to do
		elseif (isvector(adcp.(field)))
			if (length(adcp.(field)) == nf)
				% fields of length nf are not sorted,
				% because ids of files are not changed
				% adcp.(field) = adcp.(field)(sdx);
			elseif (length(adcp.(field)) == nt)
				adcp.(field)(nid) = adcp.(field);
			else
				%warning('vector of unknown length');
				warning(['vector of unknown length "', field,'"']);
				% nothing to do
			end
		elseif (2 == ndims(adcp.(field)))
			if (size(adcp.(field),1) == nf)
				% fields of length nf are not sorted,
				% because ids of files are not changed
				% adcp.(field) = adcp.(field)(sdx,:);
			elseif (size(adcp.(field),1) == nt)
				adcp.(field)(nid,:) = adcp.(field);
			else
				warning(['matrix of unknown length "', field,'"']);
				% nothing to do
			end
		elseif (3 == ndims(adcp.(field)))
			if (size(adcp.(field),2) == nt)
				adcp.(field)(:,nid,:) = adcp.(field);
			else
				% nothing to do
				warning('cube of unknown length');
			end
		else
				warning(['not processing high dimensional field "', field,'"']);
		end
	end % for field
	% must follow sorting
	adcp.fileindex  = sdx;
end % ADCP::sort

