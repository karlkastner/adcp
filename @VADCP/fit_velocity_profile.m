% Thu Jul 17 13:16:00 WIB 2014
% Karl Kastner, Berlin
%% fit velocity profile to the streamwise velocity
function obj  = fit_velocity_profile(obj,field,ensmask)
	if (~isempty(obj.velocity_profile))

	if (nargin()<3)
		ensmask = true(obj.ens.n,1);
	end

	U_bin    = obj.velocity.(field)(:,:,1);
	S_bin    = obj.bin.S;
	H        = obj.ens.H;
	binmask  = obj.mask;
	
%	if (nargin()>1)
%		velocity_profile = feval(profile_cls);
%		if (isempty(obj.velocity_profile))
%			velocity_profile = feval(profile_cls);
%		end
%	end
	
	% TODO this is not how it ought to be,
	% ensmask should be passed to vel-profile
	%param = obj.velocity_profile.param;
	
	%obj.velocity_profile.fit(U_bin(:,ensmask),S_bin(:,ensmask),H(ensmask),binmask(:,ensmask));
	%if (~isempty(obj.velocity_profile))
		obj.velocity_profile.fit(U_bin,S_bin,H,binmask,ensmask);
	%end

%	if (isempty(param))
%		param = NaN(size(obj.velocity_profile.param,1,length(ensmask)));
%	end
%	param(:,ensmask) = obj.velocity_profile.param;
%	obj.velocity_profile.param = param;

%	obj.velocity_profile.fit(U_bin,S_bin,mask,H);
%	ens.param = profile.param;
	end
end

