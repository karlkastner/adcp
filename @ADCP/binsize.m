% Mi 10. Sep 11:36:29 CEST 2014
% Karl Kastner, Berlin
%
%% bin size (vertical distance between two bins)
function binsize = binsize(obj,fid)
	if (nargin()<2)
		binsize         = ADCP.DEPTH_SLOPE*single(obj.dat.binsize);
	else
		binsize         = ADCP.DEPTH_SLOPE*single(obj.dat.binsize(fid));
	end
end % binsize

