% Mon  4 Jun 14:28:18 CEST 2018
%% convert backscatter to suspended sediment concentration
%%
%% this is the methog called "implicit" by hanes, though it is here still
%% implemented in an explicit way, as "explicit/imlicit" in hanes only
%% mean euler forward or trapezoidal integration
function [ssc, lbs_]  = backscatter2ssc_implicit(obj,R,lbs,param)
	if (nargin()<4)
		param = obj.param;
	end
	ssc = obj.ssc0(param);
	as0 = obj.as0(param);
	
	iks2 = param(1);
	xi   = param(2);

	dR = diff([0;cvec(R)]);
	nbin = size(lbs,1);

	% TODO account for attenuation along first half of cell
	% integral of ssc along distance from transducer	
	int_c_dr  = 0;
	lbs_ = NaN(size(lbs));
	for idx=1:nbin
		% compensate backscatter for attenuation by sediment
		lbs_(idx,:) = lbs(idx,:,:) + 2*xi*int_c_dr + 2*as0*dr(idx);
		% offset
		%dbs(idx,:) = [2*as*int_c_dr];
		% sediment concentration in current cell
		ssc(idx,:)  = iks2*10.^(0.1*lbs_(idx,:));
		% integral of sediment concentration along path from transducer
		int_c_dr   = int_c_dr + dR(idx)*ssc(idx,:);
	end
	ssc = ssc + ssc0;
end % backscatter2ssc

