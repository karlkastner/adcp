
meta = sanggau_metadata();
if (~exist('reload','var') || reload)
	cs_   = CrossSection();
	adcp_ = VADCP();
	Ro = [];
	us = [];
	for idx=1:length(meta.filename.discharge)
%	clear cs	
%	load(meta.filename.discharge{idx});
%	cs_(idx)   = cs;
%	clear vadcp
%	load(meta.filename.vadcp{idx});
		adcp_(idx) = VADCP(vadcp);
	end
	vadcp = adcp_;
	reload = 0;
end


