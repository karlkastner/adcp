% Wed 21 Oct 15:09:57 +08 2020
function xy2nts(obj,vadcp)
	for odx=1:length(obj)
	cdx = obj(odx).tdx;

	%dx_(idx) = quantile(dx__,0.5);
	%dy_(idx) = quantile(dy__,0.5);
	dir = obj(odx).dir;
	xyc = obj(odx).centre();

%	hyp = hypot(dx_(idx),dy_(idx));
%	dxy = [dx_(idx)/hyp;
%	       dy_(idx)/hyp]
%	R = -rot2(0*pi/2);
%	dxy = R * dxy;
%	dx(idx) = dxy(1);
%	dy(idx) = dxy(2);

	vadcp.xy2nts(cdx,xyc,dir);
	end
end

