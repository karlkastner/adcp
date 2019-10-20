% Mi 10. Sep 11:36:29 CEST 2014
% Karl Kastner, Berlin
%% convert raw transmit pulse length to SI units (m)
function lngthtranspulse = lngthtranspulse(obj)
	lngthtranspulse = ADCP.DEPTH_SLOPE*single(obj.dat.lngthtranspulse);
end % lngthtranspulse

