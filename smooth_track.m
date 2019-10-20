% Tue 29 May 14:10:11 CEST 2018
%
%% smooth a repeatedly navigated (circular) track to produce and idealized
%% average track
%
% TODO : solve as ode:
%	 from an arbitrary point on the track:
%        get kernel density estimates and second order derivative of them
%        1.) proceed to local maximum of density, from there
%        2.) follow along the ridge (i.e. along major axis of hessian)
%	 until circle is complete
%
function [XYc,S] = smooth_track(XY,XYc,ns)
	debug = false;
	fdx  = isfinite(XY(:,1));
	XY   = XY(fdx,:);

	% TODO tol for convergence
	ni = 5;
	ns = 80;
	nc = size(XYc,1);
	if (debug)
		figure(1e3+1);
		clf
	end
	for idx = 1:ni
		if (debug) plot(XYc(:,1),XYc(:,2)); hold on; end
		XYc_old = XYc;
		% TODO search object could be stored
		[id,d] = knnsearch(XYc,XY);
		%w  = 1./d.^2;
		%q = quantile(d,0.95);
		%fdx = d<q;
		np  = size(XY,1);
		fdx = true(np,1);
		w   = double(fdx);
		%w = 1./d.^4;
		if (0)
			Xc  = accumarray(id,XY(:,1).*w,[nc,1],@sum,NaN);
			Yc  = accumarray(id,XY(:,2).*w,[nc,1],@sum,NaN);
			w   = accumarray(id,w,[nc,1],@sum,NaN);
			XYc = [Xc./w,Yc./w];
		else
			%fun = @(id) wmedian(w(id),XY(id,:));
			wfun = @wmedian;
			%wfun = @wmean;
			fun = @(id) wfun(w(id),XY(id,1));
			Xc  = accumarray(id,(1:np)',[nc,1],fun,NaN);
			fun = @(id) wfun(w(id),XY(id,2));
			Yc  = accumarray(id,(1:np)',[nc,1],fun,NaN);
			XYc = [Xc,Yc];
		end
		for jdx=2:nc
			if (isnan(XYc(jdx,1)))
				XYc(jdx,:) = XYc(jdx-1,:);
			end
		end

		XYc_ = csmooth(XYc,3,[],true);
		if (idx< ni-1)
			XYc = csmooth(XYc,ns,[],true);
		else
			XYc = csmooth(XYc,5,[],true);
			%XYc = medfilt1(XYc,5);
		end
		rms(XYc_old-XYc)
	end
	if (debug) plot(XYc(:,1),XYc(:,1)); end
	% complete the track
	XYc(end+1,:) = XYc(1,:);
	if (nargout>1)
		S = cumsum([0;hypot(diff(XYc(:,1)),diff(XYc(:,2)))]);
	end
end % smooth_track

