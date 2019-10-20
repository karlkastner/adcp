% Mi 10. Sep 11:36:29 CEST 2014
% Karl Kastner, Berlin
%
%% convert the raw bin properties to si-units
function [lngthtranspulse, distmidbin1, binsize, nbins, blnk] = ...
					convert_raw_binprops_STATIC( ...
			lngthtranspulse_raw, distmidbin1_raw, binsize_raw, nbins_raw, blnk_raw)
	lngthtranspulse = ADCP.DEPTH_SLOPE*single(lngthtranspulse_raw);
	distmidbin1     = ADCP.DEPTH_SLOPE*single(distmidbin1_raw);
	binsize         = ADCP.DEPTH_SLOPE*single(binsize_raw);
	nbins           = single(nbins_raw);
	blnk		= ADCP.DEPTH_SLOPE*single(blnk_raw);
end

