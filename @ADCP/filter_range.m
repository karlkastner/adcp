% Thu Jan  8 16:06:48 CET 2015
% Karl Kastner, Berlin
%
%% filter HADCP velocity by detecting the last valid bin
%% if the bacscatter does not decreas over 10 bins, than obtstacle or intersection
function [rmask, obj] = filter_range(obj)
	nbin = size(obj.dat.VEL,1);
	emax = max(single(obj.dat.ECHO),3);
	% TODO no magic numbers
	% distance over which the intensity is expected to fall
	dn = 10;
	de = emax(dn:end,:) - emax(1:end-dn+1,:);
	rmask = ones(size(obj.dat.VEL(:,:,1)),'single');
	for idx=1:size(emax,2)
		% first cell where backscatter is increasing over the range d
		fdx = find(de(:,idx) >= 0,1);
		if (isempty(fdx))
			fdx = nbin;
		end
		% the beam is blocked at d/2
		fdx = min(nbin,fdx+dn/2);
		%lim(idx) = fdx;
		% TODO this is a test for the impact of the hadcp range
		% fdx = min(fdx,40);
		% invalidate invalid velocity bins
		rmask(fdx+1:end,idx) = NaN;
	end
%	if (apply)
%	end
	obj.rangemask = rmask;
end

