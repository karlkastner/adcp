% Sun Jul 13 11:01:30 WIB 2014
% Karl Kastner, Berlin
%% correct for platform (boat) velocity, this is the negative bed velocity
function vel = correct_for_platform_velocity_STATIC(corstr,vel,btvel,FileNumber,iflag)
%	iflag = ~iflag;
	% subtract bottom, e.g add ship velocity (v_bottom = -v_ship)
	if (nargin() < 5 || isempty(iflag) || ~iflag)
		disp('subtracting platform velocity')
		sig = -1;
	else
		disp('adding platform velocity')
		sig = 1;
	end

	field_C = {'beam','inst','ship','earth'};	
	for jdx=1:length(field_C)
		field = field_C{jdx};
		if ( isfield(vel,field) ...
			&& isfield(btvel,field) ...
			&& ~isempty(vel.(field)))
			for idx=1:4
				vel.(field)(:,:,idx) = bsxfun(@plus,vel.(field)(:,:,idx),sig*btvel.(field)(:,idx).');
			end % for idx
		end % if
	end % for jdx
end % correct_for_platfom_velocity

