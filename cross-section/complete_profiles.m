% Sun  5 Aug 09:49:05 CEST 2018
%% fill gaps in profiles
%% assumes profile to be constant in time, this is not true
%% for tidal flow in compound cross sections and near banks
function [uf,f,ubar,uu,S] = complete_profiles(u,L)
	% get pairwise profile scales
	S = zeros(n,n);

	% for each profile
	for idx=1:n(2)
	 % for succeeding profiles
	 for jdx=idx+1:n(2)
		fdx = isfinite(u(:,idx)) & isfinite(u(:,jdx));
		% u(:,i) s = u(:,jdx)
		S(idx,jdx) = u(fdx,idx) \ u(fdx,jdx);
		S(jdx,idx) = 1./S(idx,jdx);
	 end % for jdx
	end % for idx

	% complete profiles
	uu = u;
	% for each profile
	for idx=1:n(2)
		% for each position along profile
		for jdx=1:n(1)
			% fill invalid samples
			if (isnan(u(jdx,idx))
				fdx = ~isnan(u(jdx,:));
				uu(jdx,idx) = mean(S(jdx,fdx).*u(jdx,fdx));
			end % if
		end % for jdx
	end % for
	% extrapolate exponentially to side with no slip condition
	id = find(isfinite(uu(:,1)),1,'first');
	if (id > 1)
		% u = u0 (1-exp(-id/L))
		u0 = uu(id,:)/(1-exp(1/2-id/L));
		uu(1:id-1,:) = bsxfun(@times,(1-exp(1/2-(1:id-1)'/L)),u0);
	end
	id = find(isfinite(uu(:,1)),1,'last');
	if (id < n(1))
		% u = u0 (1-exp(-id/L))
		u0 = uu(id,:)/(1-exp(-(n(1)+1/2-id)/L));
		uu(id+1:n(1),:) = bsxfun(@times,(1-exp(n(1)+1/2-(id+1:n(1)))'/L)),u0);
	end
	% the velocity profile
	uf    = mean(uu,2);

	% normalized
	f    = uf/sum(uf);

	% the cross section averaged velocity
	ubar = mean(uu,1);
end % complete_profiles

