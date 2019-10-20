% Tue Oct 28 09:18:57 CET 2014
% Karl Kastner, Berlin
%% extrapolate velocity to banks, surface and bottom
%% TODO, this is only applicable for Grid2
function [U, id_] = extrapolate_velocity_2d_STATIC(U,bval,cX2,cXX,cYY,mask,clip)

	r = rmse(U,3);
	q=quantile(flat(r),[0.25 0.5 0.75]);
	% hline([q(2)+2*(q(3)-q(2))]);
	fdx = r > q(2)+2*(q(3)-q(2));
	[id, jd] = ind2sub(size(U),find(fdx));
	U(id,jd,:) = NaN;

	% do not use size here, returns only a 1x2 vector in case n3 == 1
	n1 = size(U,1);
	n2 = size(U,2);
	n3 = size(U,3);

	% zero values below bottom
	% no: this is not necessary and leads to 0 instead of nan later
%	cmsk    = bsxfun(@lt,cYY,cvec(bval));
%	cmsk    = repmat(cmsk,[1 1 n3]);
%	U(cmsk) = 0;


	% zero-valued ghost cells left and right (convex hull)
	if (1)
		cXX = [2*cXX(1,:)-cXX(2,:); cXX; 2*cXX(end,:)-cXX(end-1,:)];
		cYY = [2*cYY(1,:)-cYY(2,:); cYY; 2*cYY(end,:)-cYY(end-1,:)];
		U_  = NaN(n1+2,n2,n3);
		U_(2:end-1,:,:) = U;
		U_(1,:,:) = 0;
		U_(end,:,:) = 0;
		mask = [true(1,n2);mask;true(1,n2)];
		n1_ = n1+2;
		n2_ = n2;
	else
		cXX = [2*cXX(:,1)-cXX(:,2), cXX, 2*cXX(:,end)-cXX(:,end-1)];
		cYY = [2*cYY(:,1)-cYY(:,2), cYY, 2*cYY(:,end)-cYY(:,end-1)];
		U_  = NaN(n1,n2+2,n3);
		U_(:,2:end-1,:) = U;
		U_(:,1,:) = 0;
		U_(:,end,:) = 0;
		mask = [true(n1,1),mask,true(n1,1)];
		n1_ = n1;
		n2_ = n2+2;
	end

	% valid cells
	gdx =  isfinite(U_);
	% invalid cells above the bottom
	fdx =  bsxfun(@and,(~gdx),mask);
	% interpolate
	if (1 == n3)
		% 2d-interpolation (no time)
		% TODO scale axes
		interp  = TriScatteredInterp(cXX(gdx),cYY(gdx),U_(gdx),'linear');
		U_(fdx) = interp(cXX(fdx),cYY(fdx));
	else
		% 3d-interpolation (with time)
		% Note: this is only consistent if coefficients are values
		cXX = repmat(cXX,[1,1,n3]);
		cYY = repmat(cYY,[1,1,n3]);
		cTT(1,1,:) = 1:n3;
		cTT = repmat(cTT,[n1_, n2_, 1]);
%		cTT(1,1,1) = cTT(1,1,1)-1e-3;
%		cTT(end,1,1) = cTT(1,1,1)-1e-3;
		% nearest has to be chosen, as natural/linear does interpolate zero to the top layer
		interp  = TriScatteredInterp(cXX(gdx),cYY(gdx),cTT(gdx),U_(gdx),'nearest');
%		interp  = TriScatteredInterp(cXX(gdx),cYY(gdx),cTT(gdx),U_(gdx),'linear');
%		interp  = TriScatteredInterp(cXX(gdx),cYY(gdx),cTT(gdx),U_(gdx),'natural');
		U_(fdx) = interp(cXX(fdx),cYY(fdx),cTT(fdx));
	end
	% TODO, this is a quick fix, but linear/natural does not do a proper job
%	gdx =  isfinite(U_);
%	fdx =  ~gdx;
%	interp  = TriScatteredInterp(cXX(gdx),cYY(gdx),cTT(gdx),U_(gdx),'nearest');
%	U_(fdx) = interp(cXX(fdx),cYY(fdx),cTT(fdx));
	
	% remove convex hull
	if (1)
		U = U_(2:end-1,:,:);
	else
		U = U_(:,2:end-1,:);
	end
	id_ = [];

	if (0)

	% do not invalidate topmost cell
	mask(:,end) = true;
	mask = repmat(~mask,[1 1 n3]);
	U(mask) = NaN;

	% TODO, quick fix for top layer
	% TODO use regionfill istead of interp
	U(end,:,:) = U(end-1,:,:);
	for idx=1:size(U,2)
	 for jdx=1:size(U,3)
		if (isnan(U(end,idx,jdx)))
			U(end,idx,jdx) = 0;
		end
	 end
	end

	end
		
end % extrapolate_velocity

