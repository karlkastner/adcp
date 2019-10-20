% Mon Jul 14 00:02:55 WIB 2014
% Karl Kastner, Berlin
%
%% depth of instrument (for submerged deployments)
%
function idepth_m = instrument_depth_m(obj)
	% convert bar to instrument depth
	idepth_m = ADCP.IDEPTH_SLOPE_M_PER_BAR*obj.pressure_bar + ADCP.IDEPTH_OFFSET_M;
end

