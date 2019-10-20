% Do 5. Nov 15:10:19 CET 2015
% Karl Kastner, Berlin
%% velocity magnitude
function obj = to_abs(obj,field)
	if (nargin()<2)
		field = 'earth';
	end
	vel        = hypot( obj.velocity.(field)(:,:,1), ...
			    obj.velocity.(field)(:,:,2));
	vel(:,:,2) = atan2( obj.velocity.(field)(:,:,1), ...
			    obj.velocity.(field)(:,:,2));
	obj.velocity.abs = vel;
end

