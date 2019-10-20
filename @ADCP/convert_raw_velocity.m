% Thu Dec  4 18:59:48 CET 2014
% Karl Kastner, Berlin
%% convert scaled integer raw velocity to float SI (m/s)
function obj = convert_raw_velocity(obj)

	if (isfield(obj.dat,'btvel'))
	 	[vel_, btvel_] = obj.convert_raw_velocity_STATIC(obj.dat.VEL,obj.dat.btvel);
	else
	 	vel_ = obj.convert_raw_velocity_STATIC(obj.dat.VEL,[]);
		btvel_ = [];
		btvel  = [];
	end

	for idx=1:length(obj.dat.corstr)
		fdx = find(obj.dat.FileNumber == idx);
		switch (obj.dat.corstr{idx}(1:4))
		case {'Beam'}
			vel.beam(:,fdx,:)  = vel_(:,fdx,:);
			disp('Beam coordinates found');
		case {'Inst'}
			vel.inst(:,fdx,:)  = vel_(:,fdx,:);
			disp('Instrument coordinates found');
		case {'Ship'}
			vel.ship(:,fdx,:)  = vel_(:,fdx,:);
			disp('Ship coordinates found');
		case {'Eart'}
			vel.earth(:,fdx,:) = vel_(:,fdx,:);
			disp('Earth coordinates found');
		otherwise
			error(['Unknown coordinate system ' obj.dat.corstr{idx}]);
		end % switch
		if (~isempty(btvel_))
			switch (obj.dat.corstr{idx}(1:4))
			case {'Beam'}
				btvel.beam(fdx,:)  = btvel_(fdx,:);
			case {'Inst'}
				btvel.inst(fdx,:)  = btvel_(fdx,:);
			case {'Ship'}
				btvel.ship(fdx,:)  = btvel_(fdx,:);
			case {'Eart'}
				btvel.earth(fdx,:) = btvel_(fdx,:);
			otherwise
				error(['Unknown coordinate system ' corstr{idx}]);
			end % switch
		end
	end % for idx

	obj.velocity   = vel;
	obj.btvel = btvel;
end % convert_raw_velocity

