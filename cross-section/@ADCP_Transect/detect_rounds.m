% Tue 29 May 10:58:13 CEST 2018
%% detect rounds, i.e. when boat returns to initial position
function [d0,obj] = detect_rounds(obj,adcp)
	% must be smaller 0.5
	% if too few transects are detected, p is probably too large
	% TODO no magic number
	p = 0.3;
	x = adcp.X;
	y = adcp.Y;

	% start point
	fdx = find(isfinite(x),1);
	x0  = x(fdx);
	y0  = y(fdx);

	% diameter
	d0  = hypot(x-x0,y-y0);
	dm  = 2*nanmedian(d0);

	state = 0;
	mdx = 1;
	n = 0;
	for idx=2:length(x)
		if (~isnan(d0(idx)))
		switch (state)
		case {0} % approaching half round opposite of the beginning
			if (d0(idx) < d0(mdx))
				% new minimum
				mdx = idx;
			end
			if (d0(idx) > (1-p)*dm)
				% push back current minimum
				n = n+1;
				rid(n,1) = mdx;
				state = 1;
			end
		case {1} % approaching end of round
			if (d0(idx) < p*dm)
				mdx   = idx;
				% approach opposite end
				state = 0;
			end
		end % state
		end
	end % for idx
	% push end of last round (does not need to be completed)
	n = n+1;
	rid(n,1) = length(x);
	% mdx;
	obj.rid = rid;
end % detect_rounds

