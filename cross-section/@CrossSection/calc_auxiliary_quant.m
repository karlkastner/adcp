% Mon Aug 25 08:57:09 CEST 2014
% Karl Kastner, Berlin
%% compute auxiliary quantities
function cs = calc_auxiliary_quant(cs)
	
	% effective energy slope
	try
	u_s_     = gridN.val.u_s(:,1);
	h_       = bottom;
	slope    = nanmean(u_s_.*u_s_./h_)/Constant.g;
	cs.slope = slope;

	% effective shear velocity
	u_s      = sqrt(Constant.g*cs.radius*slope);
	cs.u_s   = u_s;

	% effective roughness length
	%ln_z0    = log(cs.depth/Constant.Euler) - Constant.KAPPA*cs.U./u_s;
	ln_z0    = log(cs.depth) - 1 - Constant.Karman*cs.U./u_s;
	cs.ln_z0 = ln_z0;

	% Chezy's coeficcient
	% van rijn: Chezy        = 5.75*sqrt(Constant.g)*log10(0.4*Rh./exp(ln_z_0));
	% sqrt(g)U/u_s == g/KAPPA ln(H/ez_0)
	Chezy    = sqrt(Constant.g)*cs.U/u_s;
	cs.Chezy = Chezy;
	% friction parapeter
	C_f      = Constant.g./Chezy.^2;
	cs.C_f   = C_f;

	catch 	e
		disp(e);
	end

end % calc_axiliary_quant

