% Sat 21 Jul 12:32:06 CEST 2018
%% convert backscatter to suspended sediment concentration
%% c.f. sassi
function ssc = backscatter2ssc_sassi(obj,R,bs,Ibs,R0,last,varargin)
	if (any(R(2,:) <= R(1,:)))
		error('R must be ascending');
	end
	% integrate the backscatter
	if (isempty(Ibs))
		nens = size(R,2);
		Ibs  = cumintL(bs,[zeros(1,nens);R]);
	end
	% backscatter at reference distance
	%bs_ref  = interp1(R,bs,R_ref,'linear');
	%Ibs_ref = interp1(R,Ibs,R_ref,'linear');
	for idx = 1:size(R,2)
		bs_ref(idx,1)  = interp1(R((1:last(idx)),idx),bs((1:last(idx)),idx),obj.R_ref,'linear','extrap');
		Ibs_ref(idx,1) = interp1(R((1:last(idx)),idx),Ibs((1:last(idx)),idx),obj.R_ref,'linear','extrap');
	end
		
	if (~isempty(R0))
		for idx = 1:size(R,2)
			bs0(idx,1)  = interp1(R((1:last(idx)),idx),bs((1:last(idx)),idx),R0(idx),'linear','extrap');
			Ibs0(idx,1) = interp1(R((1:last(idx)),idx),Ibs((1:last(idx)),idx),R0(idx),'linear','extrap');
		end
	end

	ssc = obj.backscatter2ssc_sassi_sample(obj,bs0,Ibs0,bs_ref,Ibs_ref,varargin{:});
end % backscatter2ssc_sassi

