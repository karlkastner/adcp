% Sun Jul 13 19:50:36 WIB 2014
% Karl Kastner, Berlin
%% interpolate invalid bin-samples from last and next ensemble
% TODO, this is quick and dirty
function vel = fixnan(vel)
	nens = size(vel,2);
	nbin = size(vel,1);
	for idx=1:4
		vi = vel(:,:,idx);
		fdx = find(isnan(vi));
		% exclude first and last bin
		gdx = find(fdx > nens & fdx <= nbin*(nens-1));
		vel(fdx(gdx)) = 0.5*(vel(fdx(gdx)-nens)+vel(fdx(gdx)+nens));
	end
end

