% Sun Jan  5 17:32:13 WIB 2014
% Karl Kastner, Berlin
%
%% rotate the velocity in the horizontal plane with respect to the directional
%% vector dir
%% dir : direction of the transect
function [vel] = rotate_velocity(dir,vel,iflag)
	if (isvector(dir))
		dir = cvec(dir);
	end

	% normalize
	if (1 == size(dir,2))
		dir = dir/norm(dir);
		c   = dir(1);
		s   = dir(2);
		% changed on Fri 14 Sep 14:12:38 CEST 2018
		R   = [ c, -s;
       	                s,  c];
		if (nargin()>2 && iflag)
			R = R';
		end
	else
		dir = bsxfun(@times, dir, 1./hypot(dir(1,:),dir(2,:)));
		c   = dir(1,:);
		s   = dir(2,:);
		if (nargin()>2 && iflag)
			% identical to transpose
			s = -s;
		end
	end
	siz = size(vel);
	if (3 == ndims(vel))
		U  = squeeze(vel(:,:,1));
		V  = squeeze(vel(:,:,2));
		if (isvector(dir))
			UV = [U(:), V(:)]*R.';
			vel(:,:,1) = reshape(UV(:,1),siz(1:2));
			vel(:,:,2) = reshape(UV(:,2),siz(1:2));
		else
			vel(:,:,1) =  bsxfun(@times,c,U) + bsxfun(@times,-s,V);
			vel(:,:,2) =  bsxfun(@times,s,U) + bsxfun(@times, c,V);
		end
	else
		% same as vel(:,1:2) = (R*vel(:,1:2)')';
		if (isvector(dir))
			vel(:,1:2) = (vel(:,1:2)*R.');
		else
			U = vel(:,1);
			V = vel(:,2);
			vel(:,1) =  c.*U + -s.*V;
			vel(:,2) =  s.*U +  c.*V;
		end
	end
end % rotvel

