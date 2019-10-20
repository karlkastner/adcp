% Wed Feb  5 09:16:42 WIB 2014
% Karl Kastner, Berlin
%% determine boat velocity from bottom track, inverse of bottom track

function obj = btvel_from_position(obj)
	% TODO this should be in a separate function
	dt   = Constant.SECONDS_PER_DAY*cdiff(obj.time);
	dx   = cdiff(obj.X);
	dy   = cdiff(obj.Y);

	% btvel is the negative of the ship velocity
	% TODO, this should use the vtg field, if available
	obj.utm.btvel = -[dx./dt dy./dt zeros(size(dt))];
end

