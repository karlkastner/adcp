% Mi 11. Nov 16:30:28 CET 2015
% Karl Kastner, Berlin
%
%% extrapolate backscatter to bed, surface and banks
% TODO this is only applicable to Grid2
% TODO this is a blind copy of the velocity extrapolation
% -> interpolating to zero at the bottom is not good
% -> clip can be less restrictive
function [Bs, id_] = extrapolate_backscatter_2d_STATIC(Bs,bval,cX2,cXX,cYY,mask,clip)
%	if (~isempty(clip))
%		fdx = abs(cXX(:,1) - cXX(1,1)) < clip(1);
%		Bs(fdx,:,:) = NaN;
%		fdx = abs(cXX(:,1) - cXX(end,1)) < clip(2);
%		Bs(fdx,:,:) = NaN;
%	end

	r = rmse(Bs,3);
	%q=quantile(flat(r),[0.25 0.5 0.75]);
	q=quantile(flat(r),[0.6 0.8]);
	% hline([q(2)+2*(q(3)-q(2))]);
	% accept overshoot by 100%
	fdx = r > q(2)+11*(q(2)-q(1));
	[id, jd] = ind2sub(size(Bs),find(fdx));
	Bs(id,jd,:) = NaN;

	% create a convex hull of zeros
	% bottom
	for idx=1:length(cX2)
		fdx = find(cXX(idx,:) > bval(idx),1);
		Bs(idx,fdx,:) = 0;
	end
	n1 = size(Bs,1);
	n2 = size(Bs,2);
	n3 = size(Bs,3);
	% left and right
	Bs_ = zeros(n1+2,n2,n3);
	Bs_(2:end-1,:,:) = Bs;
	cXX = [2*cXX(1,:)-cXX(2,:); cXX; 2*cXX(end,:)-cXX(end-1,:)];
	cYY = [2*cYY(1,:)-cYY(2,:); cYY; 2*cYY(end,:)-cYY(end-1,:)];

	% determine invalid cells
	gdx =  isfinite(Bs_);
	fdx =  ~gdx;
	% interpolate
	cXX = repmat(cXX,[1,1,n3]);
	cYY = repmat(cYY,[1,1,n3]);
	cTT(1,1,:) = 1:n3;
	cTT = repmat(cTT,[n1+2, n2, 1]);
%	cTT(1,1,1) = cTT(1,1,1)-1e-3;
%	cTT(end,1,1) = cTT(1,1,1)-1e-3;
	interp  = TriScatteredInterp(cXX(gdx),cYY(gdx),cTT(gdx),Bs_(gdx),'linear');
	Bs_(fdx) = interp(cXX(fdx),cYY(fdx),cTT(fdx));
	% TODO, this is a quick fix, but linear/natural does not do a proper job
	gdx =  isfinite(Bs_);
%	fdx =  ~gdx;
%	interp  = TriScatteredInterp(cXX(gdx),cYY(gdx),cTT(gdx),Bs_(gdx),'nearest');
%	Bs_(fdx) = interp(cXX(fdx),cYY(fdx),cTT(fdx));
	
	% extract
	Bs = Bs_(2:end-1,:,:);
	% do not invalidate topmost cell
	mask(:,end) = true;
	mask = repmat(~mask,[1 1 n3]);
	Bs(mask) = NaN;
	id_ = [];

end % extrapolate_backscatter_2d_STATIC

