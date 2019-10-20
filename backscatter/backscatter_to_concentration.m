% Tue 16 May 14:08:55 CEST 2017
%% convert acoustic backscatter to suspended sediment mass concentration
%% backscatter S has to be corrected for attenuation
function Cm = backstatter_to_concentration(S,ks)
	%ks2 = backscatter_coefficient(d,f);
	% average backscatter coefficient of grain diameter
	% ks2 = mean(ks2);
	Cm  = 1./ks2.*10.^(S/10);
end

