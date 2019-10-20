% Sat Jan 25 13:28:52 WIB 2014
%
%% fit the optimal cross section as the main axis of the transect
%% by regressing a line through the measurement points in the x-y plane
%%
%% y = c0 + c1 x
%
function obj = fit_cross_section(obj)

	% assure that N coordinate will be between 0 and 1 and not -1 and 0
%	if (obj.xpmax < obj.xpmin)
%		help = obj.xpmin;
%		obj.xpmin = obj.xpmax;
%		obj.xpmax = help;
%		help = obj.ypmin;
%		obj.ypmin = obj.ypmax;
%		obj.ypmax = help;
%	end

	% if transect is user defined, determine slope and intersect from that points
	% TODO numbers as enum
	switch (upper(obj.cmethod))
	case{'USER'} % cross section end points provided by the user
		%dx  = obj.xlim(:,2) - obj.xlim(:,1);
		%dy  = obj.ylim(:,2) - obj.ylim(:,1);
		%dwidth = hypot(dx,dy);
		%dir = [dx; dy]./dwidth;
		%obj.dwidth = dwidth;
		%obj.centre(1) = 0.5*(obj.xlim(:,1) + obj.xlim(:,2));
		%obj.centre(2) = 0.5*(obj.ylim(:,1) + obj.ylim(:,2));
		if (isempty(obj.T_max))
			if (isempty(obj.T_rel))
				obj.T_rel = 0.5;
			end
			obj.T_max = obj.T_rel*obj.dwidth;
		else
			obj.T_rel = obj.T_max/obj.dwidth;
		end
	case {'AUTO'} % from average and principal axis of the ensemble coordinates
		obj.transect.fit();
	end

	% write values
%	width = [];
%	dir   = [];
%	for idx=1:length(dx)
%		dwidth(idx) = hypot(dx(idx),dy(idx));
%		dir = [dx(idx); dy(idx)]./dwidth(idx);
		%sqrt(dx*dx + dy*dy);
%		obj.cs(idx).xpc = xpc;
%		obj.cs(idx).ypc = ypc;
%		obj.cs(idx).dir  = dir;
%		obj.cs(idx).dwidth = dwidth;
%	end
%	obj.fdx = fdx;
%	obj.cc  = cc;
end % function calc_cross_section()

