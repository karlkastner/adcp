% Mi 10. Sep 11:36:29 CEST 2014
% Karl Kastner, Berlin

%% convert raw distance to first bin centre to SI
function distmidbin1 = distmidbin1(obj,fid)
	if (nargin()<2)
		distmidbin1     = ADCP.DEPTH_SLOPE*single(obj.dat.distmidbin1);
	else
		distmidbin1     = ADCP.DEPTH_SLOPE*single(obj.dat.distmidbin1(fid));
	end
end % distmidbin1

