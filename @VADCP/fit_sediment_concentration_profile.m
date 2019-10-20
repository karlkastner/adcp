% Wed 17 May 14:51:21 CEST 2017
% Karl Kastner, Berlin
%% fit_suspended_sediment_concentration_profile(obj,profile_cls,ensmask,nwin)
function obj  = fit_suspended_sediment_concentration_profile(obj,ensmask)
	if (~isempty(obj.sediment_concentration_profile))
%	if (nargin()<3)
%		nwin = 0;
%	end
	% make odd
	%nwin = 1 + 2*ceil(nwin/2);

	mask  = obj.mask;
	mask  = all(mask,3);
	H     = obj.ens.H;
	Z_bin = obj.bin.Z; 
	ifield = 'sediment_concentration';
	ssc   = obj.bin.(ifield);
%	ssc   = mean(obj.bin.backscatter,3);

	% bs as ssc proxy
%	ssc   = nanmean(obj.beta,3);

%	ssc_profile = feval(profile_cls);
%	ssc_profile.nlflag = nlflag;
%	if (isempty(obj.sediment_concentration_profile))
%		obj.sediment_concentration_profile = feval(profile_cls);
%		obj.sediment_concentration_profile.nlflag = nlflag;
%	end
%	if (~isempty(obj.sediment_concentration_profile))
		%obj.sediment_concentration_profile = feval(profile_cls);
		obj.sediment_concentration_profile.fit(ssc,Z_bin,H,mask,ensmask);
%	end

	%

%	ssc   = ssc(:,ensmask);
%	Z_bin = Z_bin(:,ensmask);
%	H     = H(ensmask);
%	mask  = mask(:,ensmask);

	% TODO pass ensmask to ssc profile
%	ssc_profile.fit(ssc,Z_bin,H,mask,nwin);
%	obj.sediment_concentration_profile.c(:,ensmask) = ssc_profile.c;
%	obj.sediment_concentration_profile.s2c(:,ensmask) = ssc_profile.s2c;
	end
end % fit_suspended_sediment_profile

