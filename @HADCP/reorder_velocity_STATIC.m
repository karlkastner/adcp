% Tue Jul 29 10:45:19 WIB 2014
% Karl Kastner, Berlin
%
%% reorder the HADCP velocity data into the first three slots, the HADCP
%% has just three beams, but the software stores data for
%% four beams, similar to the four beam VADCPs
function vel = reorder_velocity(corstr,FileNumber,vel)
	for idx=1:length(corstr)
		fdx = find(FileNumber == idx);
		switch (corstr{idx}(1:4))
		case {'Beam'}
			% somehow the ADCP stores the data of beam 3 into dimension 4
			% copy to dim 4 to dim 3 and and zero dim 4 (non existing beam)
			vel.beam(:,fdx,3) = vel.beam(:,fdx,4);
			vel.beam(:,fdx,4) = 0;
		case {'Inst'}
			% zero error velocity
%			vel.inst(:,fdx,4) = 0;
		case {'Ship'}
			% zero error velocity
%			vel.ship(:,fdx,4) = 0;
		case {'Eart'}
			% somehow the ADCP stores the error velocity data  into dimension 3
			% copy dim 3 to dim 4 and zero dim 3 (vertical) 
			vel.earth(:,fdx,4) = vel.earth(:,fdx,3);
			vel.earth(:,fdx,3) = 0;
		otherwise
			error();
		end
	end
end % reorder_velocity

