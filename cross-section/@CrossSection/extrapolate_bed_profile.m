% Tue Oct 28 09:52:52 CET 2014
% Karl Kastner, Berlin
%% extrapolate bed profile to channel banks
% One could linearly extrapolate based on the gradient of the last measured part,
% but practical experience shows this not to be reliable, as river cross sections
% have approximately the form of a bathtub
function obj = extrapolate_bed_profile(obj)
	% TODO extrapolation should be part of Grid1
	% TODO use level instead of bf offse
	% TODO, this assumes constant wide bins

	val = obj.grid_n.val.zb;
	bankfull_offset = 0;
	fdx = find(isnan(val));
	gdx = find(~isnan(val));
	val(fdx,1) = interp1([0.5; gdx; length(val)+0.5],[-bankfull_offset; val(gdx,1); -bankfull_offset],fdx);
	obj.grid_n.val.zb = val;

end % extrapolate_bottom_profile


