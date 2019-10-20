% Thu 19 Jul 14:14:52 CEST 2018
% function [ssc, cbs, ac, ac0] = bs2ssc_sample(obj,bs,dIbs,R,param)
%% convert backscatter 2 suspended sediment concentration
function [ssc0, ac0, as00, bs0, Ibs0, cbs] = backscatter2ssc_sample(obj,R,bs,R0,last,varargin)
	if (1)
		switch (obj.method)
		case {'linear'}
			ssc = obj.param(1)*bs;
			next = 2;
			if (obj.withattenuation)
				Ibs  = cumintL(bs,[zeros(1,size(R,2));R]);
				ssc = ssc + obj.param(2).*bs.*Ibs(1:end-1,:);
				next = next+1;
			end
			if (obj.withbgattenuation)
				ssc = ssc + obj.param(next).*bs.*R;
				next = next+1;
			end
			if (obj.withbgconcentration)
				ssc = ssc + obj.param(next);
			end
		case {'nonlinear'}
		[ssc,ac,as0,bs_,Ibs_] = obj.backscatter2ssc(R,bs,varargin{:});
		otherwise
			error('here');
		end
		ssc0 = zeros(size(bs,2),1);
		for idx=1:size(bs,2)
			ssc0(idx,1) = interp1(R(1:last(idx),idx),ssc(1:last(idx),idx),R0(idx),'linear','extrap');
			%bs0(idx,1)  = interp1(R(1:last(idx),idx),bs_(1:last(idx),idx),R0(idx),'linear','extrap');
			%Ibs0(idx,1) = interp1(R(1:last(idx),idx),Ibs_(1:last(idx),idx),R0(idx),'linear','extrap');
			%ac0(idx,1)  = interp1(R(1:last(idx),idx),ac(1:last(idx),idx),R0(idx),'linear','extrap');
			%as00(idx,1) = interp1(R(1:last(idx),idx),as0(1:last(idx),idx),R0(idx),'linear','extrap');
		end % for idx

	else
		% this is quicker, but does not correct Ibs:
		if (nargin() < 5)
			param = obj.param;
		end
	
		np = 1;
		if (obj.withattenuation)
			np = np+1;
			b  = param(np);
			% attenuation correction factor
			ac  = 1./(1+b*dIbs);
		else
			ac = 1;
		end
	
		ssc0 = obj.ssc0(param);
		as0  = obj.as0(param);
	
		ac0  = 10.^(0.2*as0*R);
	
		% fraction sensed with backscatter
		cbs = param(1)*( bs .* ac0 .* ac );
	
		% total backscatter
		ssc = cbs + c0;
	end
end

