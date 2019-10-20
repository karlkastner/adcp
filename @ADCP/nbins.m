% Mi 10. Sep 11:36:29 CEST 2014
% Karl Kastner, Berlin
%
%% number of bins for each file
function nbins = nbins(obj,id)
	nbins = single(obj.dat.nbins);
	if (nargin()>1)
		nbins = nbins(id);
	end
end % binsize

