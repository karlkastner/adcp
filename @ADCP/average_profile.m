% Tue  5 Jun 09:54:52 CEST 2018
%
%% average backscatter for each sample within an specific interval
% TODO this is problematic when the depth changes
function [out] = average_profile(obj,val,id,t_avg,fdx)
	for idx=1:length(id)
		if (~isempty(t_avg))
			fdx = abs(obj.time - obj.time(id(idx))) <= 0.5*t_avg;
		end
		if (islogical(fdx))
			jd = find(fdx);
		else
			jd = fdx;
		end
		if (isempty(jd))
			out(idx).z = [];
			out(idx).h = [];
			out(idx).v = [];
		else
		% depth
		h  = obj.ens.H(jd);
		% position above the bottom
		z  = obj.bin.Z(jd);
		% last valid bin
		last = obj.ens.last(jd);
		% average depth
		h0    = nanmedian(h);
		% average position of first bin
		zmax  = nanmedian(z(1,:));
		% average position of last bin
		zmin  = nanmedian(z(sub2ind(size(z),cvec(last),cvec(1:length(jd)))));
		% average number of bins
		nb    = ceil(nanmedian(last));
		% refrence depth profile
		z0    = linspace(zmin,zmax,nb).';
		if (z(1,1)<z(end,1))
			z0 = flipud(z0);
		end
		
		% interpolate individual profiles to reference profile
		v0 = NaN(nb,length(jd));
		for jdx=1:length(jd)
			z0j = cvec(z2s_rational(z(1:last(jdx)),h(jdx),h0));
			vj  = val(1:last(jdx),jd(jdx));
			fdx = isfinite(z0j) & isfinite(vj);
			if (sum(fdx) > 2)
				v0(:,jdx) = interp1(z0j(fdx),vj(fdx),z0,'linear','extrap');
			end
		end % for idx

		% average the value
		v0 = nanmedian(v0,2);

		% store
		out(idx).z = z0;
		out(idx).h = h0;
		out(idx).v = v0;
		end
	end % for idx
end

