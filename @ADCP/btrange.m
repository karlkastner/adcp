% Thu Jul 17 21:04:16 WIB 2014
%
%% convert raw btrange to vertical distance (projected distance) of the bed
%% level below the transducer, when the transducer is looking vertically down
%% this is the depth less the transducer depth
function [btrange, btn] = btrange(obj)
%		if (nargin() < 2 || isempty(beamangle_rad))
%			beamangle_rad = VADCP.DEFAULT_BEAMANGLE_RAD;
%		end

	if (isempty(obj.btrange_))
	if (isfield(obj.dat,'btrange'))

	beamangle_rad = obj.beamangle_rad;

	btrange = single(obj.dat.btrange);

	% replace missing values (indicated by 0) with NAN
	btrange(0 == btrange) = NaN;

	btrange = VADCP.DEPTH_SLOPE*btrange;

	% sort in time
	[time, sdx] = sort(obj.time);
	btrange    = btrange(sdx,:);

	% for each beam
	for idx=1:4
		% interpolate missing values
		btrange(:,idx) = fixnan(time,btrange(:,idx),obj.dt_max_gap);
	end % idx

	% undo sorting
	btrange(sdx,:) = btrange;

	obj.btrange_ = btrange;
	else
		if (obj.is_facing_upward)
			 d = obj.instrument_depth_m();
			 b = obj.beamangle_rad;
			 obj.btrange_ = 1./cos(b(1))*repmat(cvec(d),1,4);
		end
%	obj.btrange_ = [];
	end
	end

	btrange = obj.btrange_;

	if (nargout() > 1)
		btn = floor((median(btrange,2) - blnk/cos(beamangle_rad))/(binsize/cos(beamangle)));
		btn = min(max(btn,1),nbins);
	end

end % btrange

