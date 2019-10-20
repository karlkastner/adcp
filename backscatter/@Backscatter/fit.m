% Thu 19 Jul 09:58:55 CEST 2018
%
%% fit backscatter coefficients
%%
%% function [res, leverage, w, obj] = fit(obj,ssc0,R0,R,bs,last,param0)
%%
%% ssc0		- ns x 1, reference concentration
%% R0            - ns x 1, distance to sample along beam
%% bs            - ns x nbin, backscatter profile per sample
%% R             - ns x nbin, distance to bin from transducer along beam
%% last          - last : index last valid bin
%% param0        - initial value for parameters
%
function [res, leverage, w, obj] = fit(obj,ssc0,R0,R,bs,last,param0)
	nc    = length(ssc0);
	ssc0  = cvec(ssc0);

	% TODO choose automatically linear / non-linear depending on flags

	switch (obj.method)
	case {'linear'}
		bs0 = zeros(nc,1);
		for idx=1:size(bs,2)
			bs0(idx,1)  = interp1(R(1:last(idx),idx),bs(1:last(idx),idx),R0(idx),'linear','extrap');
		end % for idx

		if (obj.withattenuation)
			Ibs0 = zeros(nc,1);
			Ibs  = cumintL(bs,[zeros(1,size(R,2));R]);
			for idx=1:size(bs,2)
				Ibs0(idx,1) = interp1(R(1:last(idx),idx),Ibs(1:last(idx),idx),R0(idx),'linear','extrap');
				Ibs    = cumintL(bs,[zeros(1,size(R,2));R]);
			end
		else
			Ibs0 = [];
		end
		A = obj.regmat(bs0,Ibs0,R0);
	
		% solve	
		w = ones(size(ssc0));
	
		if (1)
		d     = diag(1./sqrt(diag(A'*A)));
		param = d*((A*d) \ ssc0);
	
		%iAtA = inv(A'*A);
		iAtA = d*inv(d*A'*A*d)*d;
		else
		[Q,R] = qr(A,0);
		d     = diag(1./(diag(R)));
		param = d*((R*d) \ (Q'*ssc0));
		iAtA  = d*inv((R*d)'*(R*d))*d;
		end
		
		if (0)
			obj.param = param;
			ssc_ = [A*param, obj.backscatter2ssc_sample(R,bs,R0,last)];
			norm(ssc_(:,1)-ssc_(:,2))
			pause
		end

		res    = A*param - ssc0;
		nlflag = false;		
	case {'nonlinear'}
		% TODO weight as flag
		bs0 = zeros(nc,1);
		for idx=1:size(bs,2)
			bs0(idx,1) = interp1(R(1:last(idx),idx),bs(1:last(idx),idx),R0(idx),'linear','extrap');
		end
		w1 = bs0./ssc0;
		w  = w1./(w1+mean(w1));
		bs2ssc = @(p) obj.backscatter2ssc_sample(R,bs,R0,last,p);
		%p0 = [1e6]';
		p0 = [1e3]';
		if (obj.withattenuation)
			p0(end+1) = 0;
		end
		if (obj.withbgconcentration)
			p0(end+1) = 0;
		end
		if (obj.withbgattenuation)
			p0(end+1) = 0;
		end
		nlflag = true;
	case {'sassi'}
		% ssc = a bs/(K + bs int bs dr)
		Ibs     = cumintL(bs,[zeros(1,size(R,2));R]);
		for idx = 1:size(R,2)
			bs_ref(idx,1)  = interp1(R((1:last(idx)),idx),bs((1:last(idx)),idx),obj.R_ref,'linear','extrap');
			Ibs_ref(idx,1) = interp1(R((1:last(idx)),idx),Ibs((1:last(idx)),idx),obj.R_ref,'linear','extrap');
		end
		for idx = 1:size(R,2)
			bs0(idx,1)  = interp1(R((1:last(idx)),idx),bs((1:last(idx)),idx),R0(idx),'linear','extrap');
			Ibs0(idx,1) = interp1(R((1:last(idx)),idx),Ibs((1:last(idx)),idx),R0(idx),'linear','extrap');
		end
		bs2ssc = @(p) obj.backscatter2ssc_sassi_sample( ...
					bs0,Ibs0,bs_ref,Ibs_ref,p);

		% initial values, no attenuation
		lp = [ones(nc,1),log(bs_ref) ] \ log(bs0./ssc0);
		p0 = [exp(lp(1)); lp(2); 0]
		
		c_ = obj.backscatter2ssc_sassi_sample( ...             
                                        bs0,Ibs0,bs_ref,Ibs_ref,p0);
		nlflag = true;
	end

	if (nlflag)
		sqrt_w = sqrt(w);
		[param,resn,res,exitflag,output,lambda,A] = lsqnonlin( @(p) sqrt_w.*(bs2ssc(p) - ssc0), p0);
		% undo weighting
		res = res./sqrt_w;
		%p = param;
		% res := (A*param - ssc0);
		% A'A := 2 hessian(@(param) sum(objective(param,c,bs0,Ibs0,R0).^2), param)
		A = full(A);
		% error covariance matrix
		w_        = length(w)*w/sum(w);
		iAtA     = inv(A'*diag(w_)*A);
	end

	np      = length(param);
	serr    = sqrt(nc/(nc-np))*sqrt((res'*(w.*res))/sum(w));
	%R2      = 1-serr^2/var(ssc0);
	res_    = ssc0-sum(w.*ssc0)/sum(w);
	wvar_   = ((res_'*(w.*res_))/(sum(w)-1));
	% wvar_ = wvar(w,ssc0)
	r2.pearson    = 1-serr^2/wvar_;
	
	if (obj.extended_statistics)
		y  = ssc0;
		yp = ssc0 + res;
		np = length(param);
		r2.spearman = coefficient_of_determination(y,yp,np,'Spearman');
		r2.kendall  = coefficient_of_determination(y,yp,np,'Kendall');
		r2.hodges_lehmann = coefficient_of_determination(y,yp,np,'Hodges-Lehmann');
		mad = median(abs(res));
		ns = length(y);
		aicc = akaike_information_criterion(rms(yp-y),ns,np);
	end
	

	C        = serr.^2*iAtA;

	if (nargout > 1)
		leverage = diag(A*iAtA*A');
	end

	obj.param   = param;
	obj.C       = C;
	obj.serr    = serr;
	obj.r2      = r2;
	obj.mad     = mad;
	obj.aicc    = aicc;
	% for non-linear regression in case of an attenuation offset
	%function res = objective(param,c,bs0,Ibs0,R)
	%	cp = obj.backscatter2ssc_sample(bs0,Ibs0,R0,param);
	%
	%	res = cp-c;
	%end % objective
end

