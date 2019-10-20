% Mi 10. Sep 11:36:29 CEST 2014
% Karl Kastner, Berlin
%% blanking range, range from transduce to centre of first bin
function blnk = blnk(obj)
	blnk  = ADCP.DEPTH_SLOPE*single(obj.dat.blnk);
end % blnk

