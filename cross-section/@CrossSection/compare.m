% Fri  3 Jun 16:12:19 CEST 2016
%% interpolate for all cross-sections the values to the same time-slot
%% for comparison
function [Qt, Q, t0, cdx, msk] = compare(cs, adcp, t0)
	% interpolate all to the tame time slots
	if (nargin() < 3 || isempty(t0))
		t0 = [];
		cdx = [];
		for idx=1:length(cs)
			%valid = cs(idx).filter_transects(adcp);
			valid = cs(idx).transect.valid;
			t0  = [t0;cvec(cs(idx).transect.time(valid))];
		%	n = cs(idx).transect.n;
			n = length(cs(idx).transect.time(valid));
			cdx = [cdx; idx*ones(n,1)];
		end
		[t0 sdx] = sort(t0);
		cdx     = cdx(sdx);
	end
	Qt = [];
%	Q  = struct();
	for idx=1:length(cs)
		[Qt(:,idx) Q(idx)] = cs(idx).Qt(t0);
	end
	msk = NaN(size(Qt));
	for idx=1:length(cdx)
		msk(idx,cdx(idx)) = 1;
	end
end

