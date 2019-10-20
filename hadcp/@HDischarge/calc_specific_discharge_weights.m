% 2015-01-08 14:52:48
% Karl Kastner, Berlin
%% calculate unite discharge weights
function [cw bin] = calc_specific_discharge_weights(obj)
	h0_A = obj.h0_A;

	% specific discharge at each calibration
	switch (obj.method)
	case ('empirical')
		% specific discharge weights as measured during the calibration campaigns
		for idx=1:obj.nc
			c(:,idx) = repmat(obj.q0(idx),obj.nbin,1)./obj.qs0(:,idx);
		end
	case ('theoretical')
		q   =      (log(h0_A) - obj.ln_z0_A - 1).*h0_A.^1.5;
		Q   =      (log(obj.csh0_A) - obj.csln_z0_A - 1) .*obj.csh0_A.^1.5;
		for idx=1:obj.nc
			% theoretic discharge weights	
			% assume constant friction slope over cs for the time being
			fdx = obj.csh0_A(:,idx) > 0;
			c(:,idx) = nansum(Q(fdx,idx))./q(:,idx);
		end
	otherwise
		error('');
	end

	% discharge weights as a polynomial of stage
	% 0 : constant, 1 : linear
	A       = ones(obj.nc,obj.worder+1);
	cw      = zeros(obj.worder+1,obj.nbin);
	for idx=1:obj.nbin
		hi = (h0_A(idx,:)' - obj.mh0(idx))./obj.nh0(idx);
		for pdx=1:obj.worder
			A(:,pdx+1) = hi.^pdx;
		end
		fdx = isfinite(c(idx,:));
		cw(:,idx) = A(fdx,:) \ c(idx,fdx)';
	end
	cw     = cw';
	obj.cw = cw;
end % specific_discharge_weights

