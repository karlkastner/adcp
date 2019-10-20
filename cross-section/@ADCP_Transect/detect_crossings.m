% So 14. Sep 19:33:05 CEST 2014
% Karl Kastner, Berlin

%% detect consecutive navigation of transects (channel crossings)

function obj = detect_crossings(obj,adcp)
	if (any(diff(obj.time) < 0))
		error('Ensembles are not sorted in time, sort data sets read from multiple files by time first');
	end
	switch (obj.mode)
	case {'returning'} % single cross-section
			% paranoid : data in adcp originating from different fieles
			% are necessarily sorted in time

			%[void, sdx] = sort(obj.time);
			% TODO warn if not sorted
			%N  = adcp.N; %(sdx);
			obj.detect_crossings_returning(adcp);
			%adcp.N)
			%,transect.dwidth);
			%obj.ens.transect(sdx)    = tid;
			%obj.ens.n_transect(cdx)  = tid(end);
	case {'circling'} % for bifurcation measurements
			obj.detect_crossings_circling(adcp);
			%adcp.N(:,cdx),adcp.T(:,obj.cdx),adcp.ens.tid(:,cdx));
	otherwise
			error('unimplemented transect detection mode');
	end
end % detect_crossings

