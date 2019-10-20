% 2016-05-30 10:12:46.136138714 +0200
% Karl Kastner, Berlin
%
%% separatate individual navigation of transects,
%% for cases when the boat goes in circles and crosses the branches one after
%% the other before returning to the original cross section,
%% thus the boat does not turn at the other bank to return across the same section
%% and always navigates the cross section in the same direction
function [obj] = detect_crossings_circling(obj,adcp)
	%,N,T,tid)
	cdx = obj.tdx;

	% fetch
	tid  = adcp.ens.tid(:,cdx);
	N    = adcp.N(:,cdx);
	T    = adcp.T(:,cdx);
	%adcp.N(:,cdx),adcp.T(:,obj.cdx),adcp.ens.tid(:,cdx));

	% mark samples in transect
	% TODO hysteresis
	% TODO no magic numbers
	% TODO this is problematic if boat leaves range and returns to complete transects
	inrange = (tid ~= 0);

if (0)
	% first and last indices of transects
	first = find( 1 == diff([0; inrange]));
	last  = find(-1 == diff([inrange; 0]));
end
%	ndx  = 0;
%	next = [];
	first = [];
	last  = [];
	tdx   = 0;
	state = 0;
	C     = 0;
	% TODO, no magic number
	nhyst = 0.1*obj.dwidth;
	for idx=1:length(N)
		switch (state)
			case {0} % waiting for transect
				if (inrange(idx))
					mdx = idx;
					if (N(idx) < C)
						state = 1;
					else
						state = 2;
					end
				end
			case {1} % inside left half, last state was 0
				if (~inrange(idx))
					% partial transect, not crossing the centre
					tdx        = tdx+1;
					first(tdx) = mdx;
					last(tdx)  = idx; % should be actually max
					state      = 0;
				else
					% update, if closer towards the bank
					if (N(idx) < N(mdx))
						mdx = idx;
					end
					% boat crossed over to right half
					if (N(idx) > C+nhyst)
						tdx        = tdx + 1;
						first(tdx) = mdx;
						state     = 4;
					end
				end
			case {2} % insight right half, last state was zero
				if (~inrange(idx))
					% partial transect, not crossing the centre
					tdx        = tdx+1;
					first(tdx) = mdx;
					last(tdx)  = idx;
					state     = 0;
				else
					% update, if closer towards the bank
					if (N(idx) > N(mdx))
						mdx = idx;
					end
					% boat crosses over to left half
					if (N(idx) < C-nhyst)
						tdx = tdx+1;
						first(tdx) = mdx;
						state = 3;
					end
				end % else of if inrange
			case {3} % inside left half, last state was not 0
				if (~inrange(idx))
					% boat leaves
					last(tdx) = mdx;
					state     = 0;
				else
					% update, if closer towards the bank
					if (N(idx) < N(mdx))
						mdx = idx;
					end
					% boat crosses back to left half
					if (N(idx) > C+nhyst)
						last(tdx)  = mdx;
						tdx        = tdx+1;
						first(tdx) = mdx+1;
						state      = 4;
					end
				end
			case {4} % insight right half, last state was not 0
				if (~inrange(idx))
					% boat leaves
					last(tdx) = mdx;
					state     = 0;
				else
					% update, if closer towards the bank
					if (N(idx) > N(mdx))
						mdx = idx;
					end
					if (N(idx) < C-nhyst)
						% boat crosses back to right half
						last(tdx)  = mdx;
						tdx        = tdx+1;
						first(tdx) = mdx+1;
						state = 3;
					end
				end % else of if inrange
		end % switch state
	end % for idx
	% put last state
	switch (state)
		case {0}
			% nothing to do
		case {1,2}
			% partial transect
			tdx        = tdx+1;
			first(tdx) = mdx;
			last(tdx)  = idx;
		case {3,4}
			% close transect
			last(tdx)  = idx;
	end
%	if (state > 0)
%		ndx = ndx+1;
%		next(ndx) = mdx;
%		ndx = ndx+1;
%		next(ndx) = length(N);
%	end
%	first = cvec(next(1:2:end));
%	last  = cvec(next(2:2:end));

	% continue until first in range found
	% conitnue until half is passed, then set min
	% continue until end of range, then set max

	% TODO, quick fix
	%last = min(last,length(N));

	% write back
	obj.first = first;
	obj.last  = last;
%	obj.n     = length(first);
%	state    = 0;
%	transect = [];
%	nt = 0;
%	for idx=1:length(X)
%		switch (state)
%		case {0} % out of range
%			% range was entered
%			if (obj.inrange(X(idx),Y(idx)))
%				nt = nt+1;
%				transect(nt,1) = idx;
%				state = 1;		
%			end
%		case {1} % in range
%			% range was left
%			if (~obj.inrange(X(idx),Y(idx)))
%				transect(nt,2) = idx;
%				state = 0;
%			end
%		end % switch
%	end
%	% last sample
%	if (1 == state)
%		transect(nt,2) = length(X);
%	end
%	obj.transect_A = transect;
end % detect_crossings_at_bifurcations

