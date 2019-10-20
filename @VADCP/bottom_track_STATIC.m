% So 14. Sep 14:43:13 CEST 2014
% Karl Kastner, Berlin
%% compute bottom track coordinates
% input:
% time  : [nx1] relative time in days
% btvel : [nx2] btvel velocity in arbitraty units per second
function [X, Y, S, dX, dY, dS] = bottom_track_STATIC(time, btvel)

	% only use valid samples (TODO interpolation might yield better results)
	fdx = find(isfinite(btvel(:,1).*btvel(:,2)));
	if (~isempty(fdx))
	% time steps [s]
	dt(fdx,1) = 24*3600*[diff(time(fdx)); 0];
	% instrument movement [m]
	% the instrument moves into the opposit direction it sees the bootom moving
	dX(fdx) = -dt(fdx).*btvel(fdx,1);
	dY(fdx) = -dt(fdx).*btvel(fdx,2);
	% moved distance during one time step
	dS(fdx,1) =  sqrt(dX(fdx).*dX(fdx) + dY(fdx).*dY(fdx));
	% position relative to start point integrated from bottom velocities [m]
	X(fdx,1) = cumsum(dX(fdx));
	Y(fdx,1) = cumsum(dY(fdx));
	% travelled distance NB: slightly larger than the actual river width
	S(fdx)     = 0.5*sum(dS(fdx));

	% normalise
	X(fdx) = X(fdx) - mean(X(fdx));
	Y(fdx) = Y(fdx) - mean(Y(fdx));
	end

	fdx = find(~isfinite(btvel(:,1).*btvel(:,2)));
	X(fdx)  = NaN;
	Y(fdx)  = NaN;
	S(fdx)  = NaN;
	dX(fdx) = NaN;
	dY(fdx) = NaN;
	dS(fdx) = NaN;
end % bottom_track

