% Wed Feb  5 09:16:42 WIB 2014
% Karl Kastner, Berlin
%% filter the velocity data
% TODO make this a static function / split and rename this function
% TODO, this assumes no beam vel and errvel in beam 4
function obj = filter_velocity(obj,field)
	% fetch
	vel    = obj.velocity.(field);

	% velocity magnitude
	vmag = hypot3(vel(:,:,1),vel(:,:,2),vel(:,:,3));

	% filter mask invalid velocity samples
	velmask =  (vmag > abs(vel(:,:,4))) ...
                 & (vmag < obj.VMAX) ...
		 ... %& isfinite(obj.bin.X) ...
		 ... %& isfinite(obj.bin.Y) ...
		 & isfinite(vmag(:,:,1));

	% write back values
	obj.velmask = velmask;
	
	% save memory
	% obj.velocity.mag = vmag;
end % filter_bin

