% So 14. Sep 19:33:05 CEST 2014
% Karl Kastner, Berlin

%% groups the ensembles into transects,
%% one transect is defined as all ensembles recorded during the time the boat
%% moved from one bank to the other (return is defined as separate transect)

function obj = split_transect(obj,adcp)
	N = adcp.N;
	w = obj.dwidth();

	% TODO no magic numbers
	p = 1/6;

	nens = length(N);

	% hysteresis thresholds
	L = -w*p;
	R = +w*p;
	% centre
	C = 0;

	cid = zeros(nens,1);
	tdx = 1;
	mdx = 1;
	ldx = 1;
	state =  (N(1) > C);
	first(1) = 1;
	last = [];
	for idx=2:nens
		% partial transect detection, if there is a gap of more than 1/2 of the transect
		% concatenate broken part last transect and start anew
		if (abs(N(idx) - N(idx-1)) > 0.5*w)
			tdx = tdx+1;
			cid(ldx:idx-1) = tdx-1;
			ldx=idx;
			mdx=idx;
			state =  (N(idx) > 0.5*w);
		end

		% state machine
		switch (state)
			case {0} % first section (no flush)
			if (N(idx) > R)
				state = 2;
			end
			case {1}
			if (N(idx) < L)
				state = 3;
			end
			case {2} % look for left bank (small n)
				if (N(idx) > N(mdx))
					mdx = idx;
				elseif (N(idx) < L)
					cid(ldx:mdx)   = tdx;
					last(tdx)    = mdx;
					first(tdx+1) = mdx+1;
					ldx   = mdx+1;
					tdx   = tdx+1;
					state = 3;
				end
			case {3} % look for right bank (large n)
				if (N(idx) < N(mdx))
					mdx = idx;
				elseif (N(idx) > R)
					% the boat has entered the right side
					% partial transect detection
					cid(ldx:mdx)   = tdx;
					last(tdx)    = mdx;
					first(tdx+1) = mdx+1;
					ldx   = mdx+1;
					tdx   = tdx+1;
					state = 2;
				end
		end % switch
	end % for idx
	% last transect, covers at least 1-2p of the river width,
	% so it is kept as an individual transect
	cid(ldx:end) = tdx;
	last(end+1) = nens;

	obj.first = first;
	obj.last  = last;
	% TODO, this can be build automatically and does not need to be stored
	% obj.cid   = cid;
end % function split_transect_STATIC

