% Do 30. Jul 10:24:56 CEST 2015
% Karl Kastner, Berlin
%% average the velocity over depth
function obj = depth_average_velocity(obj,field,ensmask)
	if (nargin()<2)
		field='earth';
	end

	% fetch
	nens    = obj.ens.n;
	nbin    = obj.bin.n;
	vel     = obj.velocity.(field);
	Z       = obj.bin.Z;
	H       = obj.ens.H;
	last    = obj.ens.last;
	rangemask = obj.rangemask;
	if (isfield(obj.ens.discharge,field))
		q       = obj.ens.discharge.(field);
	else
		q = struct();
	end

	% allocate
	if (~isfield(obj.ens.velocity,field) || isempty(obj.ens.velocity.(field)))
		obj.ens.velocity.(field) = NaN(obj.ens.n,size(vel,3));	
	end

	if (nargin() < 3 || isempty(ensmask))
		ensmask = true(nens,1);
	end
	
	binmask = rangemask & all(isfinite(vel),3);

	% integrate the velocity over depth to obtain specific discharge
	q = obj.depth_integrate(q,Z,vel,binmask,H,ensmask);

	% depth averaged velocity
	obj.ens.velocity.(field)(ensmask,:) = bsxfun(@times, q.total(ensmask,:), 1./H(ensmask));

	%if (isempty(obj.ens.q))
%		obj.ens.q = 
%	if (  ~isfield(obj.ens.velocity,field) ...
%	    || isempty(obj.ens.velocity.(field)))
%		obj.ens.velocity.(field) = NaN(nbin,nens,4);
%	end
	
	% write back
	obj.ens.discharge.(field) = q;
	% obj.ens.velocity.(field)(:,ensmask,:) = vel;
end % function depth_average_velocity

