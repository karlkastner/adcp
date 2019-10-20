% Thu 19 Jul 09:34:31 CEST 2018
%% convert backscatter to suspended sediment concentration
%% c.f lee hanes / sassi, with linear relation for reference concentration
function [ssc,ac,as0,bs,Ibs] = backscatter2ssc(obj,R,bs,varargin)
	if (nargin()>3)
		param = varargin{1};
	else
		param = obj.param;
	end

	ssc0   = obj.ssc0(varargin{:});
	as0    = obj.as0(varargin{:});
	b      = obj.as(varargin{:});

	nens = size(bs,2);
	if (isvector(R))
		R = repmat(cvec(R),1,nens);
	end

	% compensate for background attenuation
	as0       = 10.^(0.2*as0*R);
	bs        = bs.*as0;

	% integral of the backscatter
	dR        = diff([zeros(1,nens); R]);

	% bs is sampled already along cell
	mbs       = bs;

	%mbs       = mid([zeros(1,nens); bs]);
	Ibs       = cumsum(dR.*mbs);

	% default, take reference at surface
	Ibs_ref = 0;

	% difference of backscatter integral with respect to that at the reference depth
	dIbs    = Ibs - Ibs_ref;

	% correction for attenuation by sensed sediment
	ac = 1./(1 + b*Ibs);

	%ssc = ks2i*( bs ./ ( 1 -  log(10)/5 * ks2i * xi * dIbeta ) );
	%ssc = obj.backscatter2ssc_sample(bs,dIbs,varargin{:});
	ssc = param(1)*bs.*ac + ssc0;
end

