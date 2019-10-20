% Do 5. Nov 15:05:53 CET 2015
% Karl Kastner, Berlin
%
%% transform velocity to cross section references
%% cs-velocity is here defined as the velocity orthogonal to the cs 
%
function obj = to_cs(obj,dir,msk)
	% velocity components point east (u) and north (v),
	% and the coordinate vectors point east (x) and north (y)
	% but for conveniens, we define the cs velocity to be swapped
	% u_cs : along n
	% v_cs : along s
	% was rotated: [-dir(2), dir 1]
	dir_ = [dir(2);dir(1)];

	if (isempty(msk))
		obj.velocity.cs   = ADCP.rotate_velocity(dir_, obj.velocity.earth,false);
		obj.btvel.cs      = ADCP.rotate_velocity(dir_, obj.btvel.earth,false);
	else
		if (~isfield(obj.velocity,'cs') || isempty(obj.velocity.cs))
			% allocate
			obj.velocity.cs   = NaN(obj.bin.n,obj.ens.n,4);
			obj.btvel.cs      = NaN(obj.ens.n,4);
		end % if
		obj.velocity.cs(:,msk,:) = ADCP.rotate_velocity(dir_, obj.velocity.earth(:,msk,:),false);
		obj.btvel.cs(msk,:) = ADCP.rotate_velocity(dir_, obj.btvel.earth(msk,:),false);
	end % if
end % to_cs

	% dir = dir/norm(dir);
	% the direction of the cross section is
	% vector from left to right bank : w*[c; s]
	% rotation matrix rotating unit vector [1; 0] (earth reference)
	% to be parallel to the cs is:
	% R(a) = [ c -s]
	%        [ s  c]
	% the matrix rotating the unit vector [1;0] to be orthogonal to the
	% cs is:
	% rot2(-90)*R(a)
	% [0,-1][c -s] = [-s,-c]
	% [1, 0][s  c]   [ c,-s]

%	% the direction orthogonal to the cross section, rotating 
%	% (rotated 90deg anticlockwise)
%%	% [ 0  1][c -s]=[-s  c]
%%	% [ 1  0][s  c] [ c  s]
%	% [ 0  1][c -s]=[ s  c]
%	% [ 1  0][s  c] [ c -s]
%	% the inverse producing streamwse vel is therefore
%	% [ s  c]^-1 = [ s  c]
%	% [ c -s]      [ c -s]

