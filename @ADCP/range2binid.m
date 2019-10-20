% Sat 26 May 14:01:29 CEST 2018
%% convert distance to transducer to bin index
%
function id = range2binid(obj,r)
	r0  = obj.distmidbin1(1);
	nb  = max(obj.nbins);
	dr  = obj.binsize(1);
	id  = round((r-r0)/dr);
	id  = min(max(id,1),nb);
end

