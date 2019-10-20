% Sun 28 Jan 13:10:55 CET 2018
% Karl Kastner, Berlin
%
%% depth integrate the velocity to obtain specific discharge
%
% TODO use last, not mask
% TODO set invalid values to 0 and just decrease n
% TODO fieldname should be sufficient as argument
function q = depth_integrate(q,Z,vel,binmask,H,ensmask)
	if (isempty(q))
		q = struct();
	end

	if (nargin() < 6 || isempty(ensmask))
		nens = size(vel,2);
		ensmask = true(nens,1);
	end

	% allocate memory
	n3                     = size(vel,3);
	q.top(ensmask,1:n3)    = NaN;
	q.mid(ensmask,1:n3)    = NaN;
	q.bottom(ensmask,1:n3) = NaN;
	q.total(ensmask,1:n3)  = NaN;

	edx = find(ensmask);

	% for each ensemble
	for idx=rvec(edx)
		% only use valid bins
		vi    = vel(binmask(:,idx),idx,:);
		zi    = Z(binmask(:,idx),idx);
		n     = size(vi,1);
		if (n>0)
			% linearly extrapolate top-velocity
			% does not require no-shear condition to hold
			% defaults to constant extrapolation when only one bin is valid
			id2          = min(2,n);
			vtop         = vi(1,:) + (H(idx)-zi(1))*(vi(1,:)-vi(id2,:))/(zi(1)-zi(id2));

			% trapezoidal rule for top flow
			q.top(idx,:) = 0.5*(vtop+vi(1,:)).*(H(idx)-zi(1));

			% trapezoidal rule for mid flow
			if (n>1)
				q.mid(idx,:) = 0.5*sum(bsxfun(@times,(vi(1:end-1,:)+vi(2:end,:)),-diff(zi)));
			else
				q.mid(idx,:) = 0;
			end

			% no-slip bc for bottom flow (vbot=0) (also trapezoidal rule)
			% TODO this is not true for the concentration, only for fluxes
			q.bottom(idx,:) = 0.5*vi(end,:)*zi(end);
		end
	end % for idx

	% total discharge
	q.total(ensmask,:) =      q.top(ensmask,:) ...
				+ q.mid(ensmask,:) ...
				+ q.bottom(ensmask,:);
end

