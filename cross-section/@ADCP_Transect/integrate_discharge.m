% Mon 30 May 13:48:25 CEST 2016
% Karl Kastner, Berlin
%
%% integrate discharge
%%
%% Q  = sum q
%% q  = A_n*u_s = h dn us
%%    = h * [dx, dy]*[-v; u]
%%    = h * dt * [-ub, -vb] * [-v; u]
%%
%% note that uvb * dt is usually more accurate than dx of GPS position,
%% if uvb determined by doppler shift of ADCP bottom echo,
%% except when the GPS position (or velocity) is determined from the carrier frequency
%%
%% note that projection can be left out, if cs is defined with transect individual end points,
%% but not recommended, if there are strong secondary currents as encountered at
%% bends or bifurcations
%%
function obj = integrate_discharge(obj,adcp)
	dw_max     = 2;
	pw_min     = 0.5;

	% fetch
	time      = Constant.SECONDS_PER_DAY*adcp.time;
	UV        = adcp.ens.velocity.earth(:,1:2);
	UVb       = adcp.btvel.earth(:,1:2);
	H         = adcp.ens.H;
	dir       = obj.dir;
	ensQ      = adcp.ens.discharge.earth;

%	first     = obj.first_valid;
%	last      = obj.last_valid;

	% first and last ensemble
	first = obj.first;
	% TODO quick fix
	last  = obj.last-1;

	time_     = [];
	nan = NaN(obj.n,1);
	area      = struct('centre',nan,'left',nan,'right',nan);
	discharge = struct('total',nan,'mid',nan,'top',nan,'bottom',nan,'centre',nan,'left',nan,'right',nan);
	isvalid   = [];
	velocity = NaN(obj.n,4);

	% for each transect
	for idx=1:obj.n
		% find first and last and valid ensemble
		flag = all(isfinite([UV(first(idx):last(idx),:), ...
			             H(first(idx):last(idx)), ...
				     UVb(first(idx):last(idx),1:2)]),2);

		first(idx) = first(idx) + find([flag;1],1,'first')-1;
		last(idx)  = first(idx) + find([1;flag],1,'last')-1;
		first(idx) = max(1,first(idx));
		last(idx)  = min(last(idx),adcp.ens.n);

		%for idx=1:obj.n
		fdx  = (first(idx):last(idx))';
		fdx_ = all(isfinite(UV(fdx,1:2)),2) ...
			& all(isfinite(UVb(fdx,1:2)),2);
		fdx  = fdx(fdx_);

		if (length(fdx) < 2)
			% there must be at least two valid ensembles in a transect
			%Q = NaN;
			%A = NaN;
			t = NaN;
			W_lim = 0;
			W = NaN;
			rmseQ = NaN;
		else
			% time of transect
			t = mean(time(fdx));
	
			% time step
			dt  = cvec(cdiff(time(fdx)));

			% space step
			% TODO alternatively use coordinates from GPS
			dx = bsxfun(@times,dt,UVb(fdx,:));
			
			% space step projected to cross section
			for field_ = {'total','mid','top','bottom'}
			field = field_{1};

			% discharge per ensemble
			%if (strcmp(field,'total'))
			%dq.total =  H(fdx).*(  dx(:,1).*UV(fdx,2) ...
			%	             - dx(:,2).*UV(fdx,1) );
			%else
			dq.(field) =   dx(:,1).*ensQ.(field)(fdx,2) ...
				     - dx(:,2).*ensQ.(field)(fdx,1);
			%end

			% parallel flow
			dp.(field) =   dx(:,1).*ensQ.(field)(fdx,1) ...
				     + dx(:,2).*ensQ.(field)(fdx,2);

			%dq.top =    dx(:,1).*ensQ.top(fdx,2) ...
			%	  - dx(:,2).*ensQ.top(fdx,1);

			%dq.bottom =   dx(:,1).*ensQ.bottom(fdx,2) ...
			%	    - dx(:,2).*ensQ.bottom(fdx,1);
		
			% discharge through part of cross section

			if (0) % projected discharge
			q =  dt.*H(fdx) ...
				      .* ( UV(fdx,1).*dir(2) ...
			                 - UV(fdx,2).*dir(1) ) ...
				      .* ( UVb(fdx,1).*dir(1) ...
					 + UVb(fdx,2).*dir(2) );
			end
	
			% TODO no magic numbers
			%q = medfilt1(double(q),7);
			% Trapezoidal rule
			% q = 0.5.*(q(1:end-1)+q(2:end));
		
			discharge.(field)(idx,1)  = sum( dq.(field) );
			discharge.(field)(idx,2)  = sum( dp.(field) );
			%qparallel.(field)(idx,1)  = sum( dp.(field) );

			end % for field

			% average velocity
			% velocity = int q dt / int h dt
			velocity(idx,:) = bsxfun(@times,sum(dt.*ensQ.total(fdx,:)),1./sum(dt.*H(fdx)));

			% error estimate
			d2q   = cdiff(dq.total,2);
			rmseQ = 1/2*rms(d2q.^2);

			% dw
			dw = ( dx(:,1).*dir(1) + dx(:,2).*dir(2) );
			%dw = ( dx(:,1).*dir(2) + dx(:,1).*dir(1) );

			% fraction of cross section covered
			W     = sum(dw);
			W_lim = sum(min(dw_max,max(-dw_max,dw)));

			% projected area per ensemble
			da = H(fdx) .* dw;

			area.centre(idx) = sum(da);
		end % if two valid samples

		W     = abs(W);
		W_lim = abs(W_lim);
		time_(idx)              = t;
		discharge.rmse.centre(idx,1) = rmseQ;
		width(idx)              = W;
		width_limited(idx)      = W_lim;
		isvalid(idx)            = W_lim > pw_min.*obj.dwidth;
	end % for idx

	% make sign of area positive
	sig = sign(area.centre);
	for field_ = {'total','mid','top','bottom'}
		field = field_{1};
		discharge.(field) = bsxfun(@times,discharge.(field),sig);
		%qparallel.(field) = qparallel.(field)./sig;
	end
	area.centre = abs(area.centre);

	% relabel
	discharge.centre = discharge.total;
	discharge = rmfield(discharge,'total');
	%qparallel.centre = qparallel.total;
	%qparallel = rmfield(qparallel,'total');

	% make the time vector continuous
	fdx = isnan(time_);
	idx = 1:length(time_);
	if (sum(~fdx)>1)
		time_(fdx) = interp1(idx(~fdx),time_(~fdx),idx(fdx),'linear','extrap');
	end

	% write back
	obj.first_valid   = first;
	obj.last_valid    = last;
	obj.discharge     = discharge;
	%obj.discharge.parallel = qparallel;
	obj.velocity	  = velocity;
	obj.area          = area;
	obj.isvalid       = isvalid;
	obj.width         = width;
	obj.width_limited = width_limited;
	obj.time          = time_/Constant.SECONDS_PER_DAY;
	obj.sig           = sig;

	obj.extrapolate_to_bank(adcp);
end % integrate_discharge

