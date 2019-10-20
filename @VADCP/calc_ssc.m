% Thu 19 Jul 09:17:35 CEST 2018
%% calculate the backscatter
function obj      = calc_ssc(obj)
	bs        = obj.backscatter;
	R         = obj.R;
	bs        = obj.backscatter;
	% TODO not valid for HADCPs (trhee beams only)
	bs  = mean(bs,3);
	% TODO, for each file separately
	R = R(:,1);
	% TODO filtering?
	obj.bin.sediment_concentration   = obj.backscatter_obj.backscatter2ssc(R,bs);
end

