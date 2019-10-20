% Sun 27 May 15:41:29 CEST 2018
% TODO this is superfluous
function obj = assign_crossings(obj,transect)
	% for each transect
	for tdx=1:length(transect)
		first = transect(idx).first;
		last  = transect(idx).last;

		obj.ens.n_transect(tdx) = length(first);
		% for each crossing within the transect
		for idx=1:obj.ens.n_transect(tdx)
			fdx = first(idx):last(idx);
			obj.ens.transect(tdx,fdx) = idx;
		end
	end
end % assign_crossings

