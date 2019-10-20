% Sun 14 May 13:35:43 CEST 2017
% Karl Kastner, Berlin
%% depth integrated sediment discharge
function obj = depth_integrate_sediment_discharge(obj,field,varargin)
	% fetch
	vel     = obj.velocity.(field);
	% linear backscatter
	bs      = obj.backscatter;
	Z       = obj.bin.Z;
	H       = obj.ens.H;
%	last    = obj.ens.last;
	rangemask    = obj.rangemask;
%	qbs     = obj.ens.qbs;
%	bs_integrated = obj.ens.bs_integrated;

	% average backscatter over all beams
	% TODO HADCP only 3 beams
	% bs = mean(bs,3);
	ifield = {'sediment_concentration','backscatter'};
	ofield = {'sediment_discharge','qbs'};

	for idx=1:length(ifield)

	val     = obj.bin.(ifield{idx});

	binmask = rangemask & isfinite(val);

	% integrate the ssc over depth
	% TODO cbar
	% obj.ens.bs_integrated = depth_integrate(obj.ens.bs_integrated,Z,bs,binmask,H,varargin{:});
	%obj.ens.bs_integrated = depth_integrate(obj.ens.bs_integrated,Z,bs,binmask,H,varargin{:});

	% sediment flux in all flow directions
	flux    = bsxfun(@times, val, vel);
	
	binmask = rangemask & all(isfinite(flux),3);

	% integrate the ssc flux over depth
	obj.ens.(ofield{idx}) = obj.depth_integrate(obj.ens.(ofield{idx}),Z,flux,binmask,H,varargin{:});
	
	end
end % function integrate_backscatter_flux

