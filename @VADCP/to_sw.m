% Do 5. Nov 15:10:19 CET 2015
% Karl Kastner, Berlin
% TODO also btvel and ensmask
%% transform velocity with respect to depth averaged streamwise velocity
function obj = to_sw(obj)
	obj.velocity.sw = rotate_velocity_sw(obj.velocity.earth, obj.mask);
end

