% Fri  1 Jun 12:53:52 CEST 2018
%% extrapolate value beyond end of cross section
function v = extrapolate_n(obj,N,v)
	 % for each time slice
	 for jdx=1:size(v,2)
		fdx = isfinite(v(:,jdx));
		if (sum(~fdx) > 0 && sum(fdx)>1)
			v(~fdx,jdx) = interp1(N(fdx),v(fdx,jdx),N(~fdx),'linear','extrap');
		end % if
	 end % for jdx
end % extrapolate_S

