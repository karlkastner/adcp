% Wed Nov 12 10:41:23 CET 2014
% Karl Kästner, Berlin
%
%% rotate velocity to local streamwise reference
%% input velocity can have arbitrary reference
function vel_sw = rotate_velocity_sw(vel,valid)
%	% allocate memory
%	vel_sw = vel; %zeros(size(vel),class(vel));

	% TODO 
	% this only uses the measured cells and thus slightly differs from
	% the mean velocity over the entire depth
	% use actual depth averaged velocity

	u = vel(:,:,1);
	v = vel(:,:,2);
	if (nargin() > 1 && ~isempty(valid))
		u(~valid) = NaN;
		v(~valid) = NaN;
	end

	u = nanmean(u); %vel(:,:,1).*valid);
	v = nanmean(v); %vel(:,:,2).*valid);
	vel_sw = ADCP.rotate_velocity([rvec(u); rvec(v)],vel,true); 
if (0)
	%ih = 1./sqrt(u.*u + v.*v);
	ih = 1./hypot(u,v);
	s  = repmat(v.*ih,size(vel,1),1);
	c  = repmat(u.*ih,size(vel,1),1);
	vel_sw(:,:,1) =  c.*vel(:,:,1) + s.*vel(:,:,2);
	vel_sw(:,:,2) = -s.*vel(:,:,1) + c.*vel(:,:,2);

	vel_sw_ = ADCP.rotate_velocity([rvec(u); rvec(v)],vel);
	nannorm(flat(vel_sw(:,:,1:2)-vel_sw_(:,:,1:2)))
	vel_sw_ = ADCP.rotate_velocity([rvec(u); rvec(v)],vel,true);
	nannorm(flat(vel_sw(:,:,1:2)-vel_sw_(:,:,1:2)))
pause
end


%	for idx=1:size(vel,2)
%		u = nanmean(vel(1:ldx(idx),idx,1));
%		v = nanmean(vel(1:ldx(idx),idx,2));
%		ih = 1/sqrt(u*u+v*v);
%		s = v*ih;
%		c = u*ih;
%		vel_sw(:,idx,1) =  c*vel(:,idx,1) + s*vel(:,idx,2);
%		vel_sw(:,idx,2) = -s*vel(:,idx,1) + c*vel(:,idx,2);
%	end
end % rotvel_sw

