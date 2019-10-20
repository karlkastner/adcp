% Thu Jul 17 21:26:14 WIB 2014
% Karl Kastner, Berlin
%% unprojected (slanted) distance between the transducer and cell centres
function R = R(obj)
	beamangle_rad = obj.beamangle_rad;
	% TODO, for hadcp systems, the range of the central beam is D !!!
	% distance from transducer to cell rentres along the beams
	% dr         = binsize/cos(beamangle_rad);
	Dt = obj.Dt;
	R  = zeros(size(Dt));
	for idx=1:size(Dt,2)
		R(:,idx)  = (1/cos(beamangle_rad(idx)))*Dt(:,idx);
	end
end

