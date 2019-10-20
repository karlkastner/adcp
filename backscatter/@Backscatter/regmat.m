% Thu 19 Jul 09:58:55 CEST 2018
%
%% regression matrix
%
% ssc = p1 bs/(1 - int p2 bs + p3 dr) + p4
%     ~ p1 bs + p2' bs int bs dr + p3' bs r + p4
function A = regmat(obj,bs0,Ibs0,R0)
	% construct regression matrix
	A      = bs0;
	if (obj.withattenuation)
		A(:,2)     = bs0.*Ibs0;
	end

	if (obj.withbgattenuation)
		A(:,end+1) = bs0.*R0;
	end

	if (obj.withbgconcentration)
		A(:,end+1) = 1;
	end
end % regmat

