% Fri  1 Jun 12:53:52 CEST 2018
%% extrapolate missing values along the vertical
function v = extrapolate_S(obj,N,v)
	% for each sample along the vertical
	for idx=1:size(v,1)
	 % for each time slice
	 for jdx=1:size(v,3)
		fdx = isfinite(v(idx,:,jdx));
		if (sum(~fdx) > 0 && sum(fdx) > 1)
			v(idx,~fdx,jdx) = interp1(N(fdx),v(idx,fdx,jdx),N(~fdx),'linear','extrap');
		end % if
	 end % for jdx
	end % for idx
end % extrapolate_S

