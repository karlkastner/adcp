% Tue Jul 29 08:21:57 WIB 2014
% Karl Kastner, Berlin
%% correct RDI HADCP firmware bug (2014)
%% this bug successively invalids data every 4th-bin, which led to 3-beam solutions
%% and consequentially jumps of the transformed velocities
function [vel, fixed] = firmware_fix_STATIC( corstr,FileNumber,firmver,firmrev,vel,...
			btvel,heading_rad,pitch_rad,roll_rad,beamangle_rad)
	if (nargin()<10 || isempty(beamangle_rad))
		beamangle_rad = HADCP.DEFAULT_BEAMANGLE_RAD;
	end

	% FileNumber ist not necessarily sorted
	for idx=1:max(FileNumber)
		fdx = find(FileNumber == idx);
		% check is this firmware version is affected
		if (firmver(idx) == 11 && firmrev(idx) >= 9 && firmrev(idx) <= 10)
			if (~isempty(btvel))
				btvel_i = btvel(fdx,:);
			else
				btvel_i = [];
			end

			switch ((corstr{idx}(1:4)))
			case {'Eart'}
				% fix earth vel by filtering
				X = (0:size(vel.earth,1)-1)'; 
				X = X-mean(X);
				X = X/max(X);
				for jdx=fdx
				 for ddx=1:2
					Y = vel.earth(:,jdx,ddx);
					gdx = isfinite(Y);
					A = vander_1d(X(gdx),32);
					A = orth(A);
					Y1 = A*(A\Y(gdx));
					vel.earth(gdx,jdx,ddx) = Y1;
			 	 end
				end
		
%			% convert velocities to beam velocities first
%			if (strcmp(corstr{idx}(1:4),'Eart'))
%				[vel.ship(:,fdx,:) btvel_i] = ADCP.ship_to_earth_STATIC(vel.earth(:,fdx,:),btvel_i,heading_rad(fdx),1);
%				flag = 1;
%			else
%				flag = 0;
%			end
%			if (flag | strcmp(corstr{idx}(1:4),'Ship'))
%				[vel.inst(:,fdx,:) btvel_i] = ADCP.instrument_to_ship_STATIC(vel.ship(:,fdx,:),btvel_i,pitch_rad(fdx),roll_rad(fdx),1);
%				flag = 1;
%			end
%%			vel_s(:,fdx,:) = vel(:,fdx,:);
%			if (flag | strcmp(corstr{idx}(1:4),'Inst'))
%				% save old vel
%				ovel = vel;
%				% convert to beam velocities
%				bdx = 1:4:size(vel.inst,1);
%				% individual 2 and 3 beam solutions
%				[vel.beam(bdx,fdx,:) btvel_i] = HADCP.beam_to_instrument_STATIC(ovel.inst(bdx,fdx,:),btvel_i,beamangle_rad,123,1);
%				bdx = 2:4:size(vel.inst,1);
%				[vel.beam(bdx,fdx,:) btvel_i] = HADCP.beam_to_instrument_STATIC(ovel.inst(bdx,fdx,:),btvel_i,beamangle_rad,12,1);
%				bdx = 3:4:size(vel.inst,1);
%				[vel.beam(bdx,fdx,:) btvel_i] = HADCP.beam_to_instrument_STATIC(ovel.inst(bdx,fdx,:),btvel_i,beamangle_rad,13,1);
%				bdx = 4:4:size(vel.inst,1);
%				[vel.beam(bdx,fdx,:) btvel_i] = HADCP.beam_to_instrument_STATIC(ovel.inst(bdx,fdx,:),btvel_i,beamangle_rad,23,1);
%			end
%			if (~isempty(btvel))
%				btvel(fdx,:) = btvel_i;
%			end
%%			vel_b(:,fdx,:) = vel(:,fdx,:);
%			% fix pattern (earth coordinates only)
%			if (0)
%				% this does not work well when the velocity changes abruptly
%				for jdx=1:3
%				pattern(:,fdx,jdx) = [ zeros(4,length(fdx));
%					               0.5*(vel.beam(1:end-8,fdx,jdx)+vel.beam(9:end,fdx,jdx));
%				                       zeros(4,length(fdx)) ];
%				end
%				pattern(:,fdx,4) = 0;
%				mu_(:,fdx,1)        = 1/3*(pattern(1:4:end-1,fdx,1) + pattern(2:4:end,fdx,1) + pattern(3:4:end,fdx,1));
%				mu_(:,fdx,2)        = 1/3*(pattern(1:4:end-1,fdx,2) + pattern(2:4:end,fdx,2) + pattern(4:4:end,fdx,2));
%				mu_(:,fdx,3)        = 1/3*(pattern(1:4:end-1,fdx,3) + pattern(3:4:end,fdx,3) + pattern(4:4:end,fdx,3));
%				mu_(:,fdx,4) 	   = 0;	
%				for jdx=1:4
%					mu_i = reshape(mu_(:,fdx,jdx),1,[]);
%					mu_i = repmat(mu_i,4,1);
%					mu(1:4*size(mu_,1),fdx,jdx) = reshape(mu_i,4*size(mu_,1),[]);
%				end
%				mu(149,:,end) = 0;
%				% subtract pattern and compensate for mean
%				d = pattern-mu;
%				vel.beam(:,fdx,:) = vel.beam(:,fdx,:) - pattern(:,fdx,:) + mu(:,fdx,:);
%			end
%
			% interpolate missing values
%			if (flag | strcmp(corstr{idx}(1:4),'Beam'))
			case {'Beam'}
			% central
			vel.beam(2:4:end,fdx,3) = 0.5*(vel.beam(1:4:end-1,fdx,3) + vel.beam(3:4:end,fdx,3));
			vel.beam(3:4:end,fdx,2) = 0.5*(vel.beam(2:4:end,fdx,2)   + vel.beam(4:4:end,fdx,2));
			vel.beam(4:4:end,fdx,1) = 0.5*(vel.beam(3:4:end,fdx,1)   + vel.beam(5:4:end,fdx,1));
			% group average
			%vel(4:4:end,fdx,1) = 1/3*(vel(1:4:end-1,fdx,1) + vel(2:4:end,fdx,1) + vel(3:4:end,fdx,1));
			%vel(3:4:end,fdx,2) = 1/3*(vel(1:4:end-1,fdx,2) + vel(2:4:end,fdx,2) + vel(4:4:end,fdx,2));
			%vel(3:4:end,fdx,2) = 2*vel(2:4:end,fdx,2) - vel(1:4:end-1,fdx,2);
			%vel(3:4:end,fdx,2) = 2*vel(4:4:end,fdx,2) - vel(5:4:end,fdx,2);

			% somehow is the velocity in earth coordinates even more corrupted than in 
			% beam coordinates, so the high frequent changes are filtered here
			% TODO maybe the results are better when this is done before interpolation
			% TODO take care of boundary values
			if (0)
				velo = vel.beam;
				l = -3;
				r =  4;
				x = (l:r)'/r;
				A = [x.^0 x x.^2 x.^3]; % x.^4 x.^5 x.^6]; % x.^7]; % x.^8 x.^9 x.^10];
				Ai = inv(A'*A)*A';
				for ddx=1:3
				for jdx=(1-l):size(vel,1)-r
					p = Ai *  squeeze(velo(jdx+l:jdx+r,fdx,ddx));
					vel.beam(jdx,fdx,ddx) = p(1,:);
				end
				end
			end
			fixed(idx) = true;
			otherwise
				fixed(idx)=false;
			end % switch
		else
			fixed(idx) = false;
		end % if firmver
	end % for idx
end % function

