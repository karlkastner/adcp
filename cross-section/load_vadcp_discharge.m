% Sun Nov  2 11:31:32 CET 2014
% Mon Feb  2 14:44:42 CET 2015
% Karl Kastner, Berlin
%
%% load previously computed vadcp discharge (auxiliary function for plotting)
%
%% this function stacks data from several vadcp reference measurements into one structure
%% This assumes that all data sets where processed with the same settings
function calib = load_vadcp_discharge(opt,vflag)
	if (nargin() < 2 || isempty(vflag))
		vflag = false;
	end
	calib             = struct();
	calib.filename_C  = opt.filename.vadcp;
	calib.nc          = length(calib.filename_C);
	opt_str           = CrossSection.optstr(opt);
	n                 = 0;
	calib.time        = [];
	calib.zs0         = [];
	calib.area        = [];
	calib.width       = [];
	calib.radius      = [];
	calib.perimeter   = [];
	calib.Q           = [];
	calib.ubar        = [];
	calib.vbar        = [];
	calib.U           = [];
	calib.zb.A        = [];
	calib.u_s.A.val   = [];
	calib.ln_z0.A.val = [];
	calib.perturbation.A.val = [];
	firsttime         = true;
	for idx=1:calib.nc
		if (length(calib.filename_C{idx}) > 0)
		iname     = [calib.filename_C{idx}(1:end-4),opt_str,'.mat'];
		fprintf(1,'loading %s\n',iname);
		load(iname,'cs');
				
		calib.cs_(idx,1)       = cs;
		
		% time varying quantities, cross sectionally averaged
		t0               = cs.t0;
		if (isempty(cs.dt))
			t0 = head(cs.t0);
		end
		nt               = length(t0);
		calib.id{idx}    = cvec(n+1:n+nt);
		calib.time(n+1:n+nt,1)  = t0;
		calib.tstart(idx)=t0(1);
		calib.zs0(n+1:n+nt,1)    = cs(1).level_t(t0);
		calib.area(n+1:n+nt,1)   = cs.area(t0);
		calib.width(n+1:n+nt,1)  = cs.width(t0);
		calib.radius(n+1:n+nt,1)  = cs.radius(t0);
		calib.perimeter(n+1:n+nt,1)  = cs.perimeter(t0);
		calib.Q(n+1:n+nt,1) = cs.discharge(t0);
		calib.U(n+1:n+nt,1) = cs.U_t(t0);

		calib.zb.A(:,idx)   = cs.zb; 

		% local quantities (neither averaged over time nor cross section)
		calib.U_tn(:,n+1:n+nt) = cs.U_tn(t0);
		
		calib.q_tn(:,n+1:n+nt) = cs.q_tn(t0);

		% roughness
		calib.ln_z0.A.val(:,n+1:n+nt) = cs.ln_z0_tn(t0);
%		calib.ln_z0.A.err(:,idx) = cs.gridN.err.ln_z0(:,1);

		% shear
		calib.u_s.A.val(:,n+1:n+nt)   = cs.us_tn(t0);
%		calib.u_s.A.err(:,idx)   = cs.gridN.err.u_s(:,1);

		% wake parameter
		p = cs.var_tn('perturbation',t0);
		calib.perturbation.A.val(:,n+1:n+nt)  = p; 
		%cs.perturbation_tn(t0);
		% TODO H = zs-zb !!!
		calib.wwake_tn(:,n+1:n+nt) = bsxfun(@times,p,1./(-cs.zb));

		% linear profile parameter
		%calib.vertical_profile_parameter(:,:,n+1:n+nt) = cs.gridNr.vali.param;
%		calib.vertical_profile_parameter.a(:,n+1:n+nt) = reshape(cs.gridNr2.vali.a,cs.gridNr2.n)';
%		calib.vertical_profile_parameter.b(:,n+1:n+nt) = reshape(cs.gridNr2.vali.b,cs.gridNr2.n)';
%		calib.vertical_profile_parameter.c(:,n+1:n+nt) = reshape(cs.gridNr2.vali.c,cs.gridNr2.n)';
%		calib.vertical_profile_parameter.b(:,n+1:n+nt) = cs.gridNr2.vali.b;
%		calib.vertical_profile_parameter.c(:,n+1:n+nt) = cs.gridNr2.vali.c;

		% was qbs
		calib.sediment_discharge_tn(:,n+1:n+nt)   = cs.var_tn('sediment_discharge',t0);
		
		calib.rouse_tn(:,n+1:n+nt) = cs.rouse_tn(t0);

		if (vflag)
		if (size(cs.gridNZ.val.U,3) == length(cs.itime)-1)
			% quick fix for zero order
			U = cs.gridNZ.val.U;
			% fix different depth
			n1 = size(U,2);
			n2 = size(calib.U,2);
			if (size(U,2) > size(calib.U,2))
			elseif (size(U,2) < size(calib.U,2))
				U_ = U;
				U = NaN(size(calib.U,1),size(calib.U,2),size(U_,3));
				U(:, n2-n1+1:n2, :) = U_;
			else
			end
			calib.U(:,:,n_+1:n-1) = U;
			%calib.U(:,:,end+1:n-1) = cs.gridNZ.val.U;
			calib.U(:,:,n) = NaN;
	%		calib.Uerr(:,:,n) = NaN;
		else
 			U = cs.gridNZ.val.U;
			% fix different depth                                   
                        n1 = size(U,2);                                         
                        n2 = size(calib.U,2);                                   
                        if (size(U,2) > size(calib.U,2))                        
                        elseif (size(U,2) < size(calib.U,2))                    
                                U_ = U;                                         
                                U = NaN(size(calib.U,1),size(calib.U,2),size(U_,3));
                                U(:, n2-n1+1:n2, :) = U_;                       
                        else                                                    
                        end
			calib.U(:,:,n_+1:n) = U;
			%calib.Uerr(:,:,n_+1:n) = cs.gridNZ.err.U;
		end
%		calib.Uerr(:,:,idx) = cs.gridNZ.err.U;
		end % if vflag

		%calib.m1(:,idx) = cs(1).gridN.id.i1.m;

		if (firsttime)
			% averaged and constant quantities, taken from first file
			% assumed constant (e.g. same processing setting for all data sets)
			calib.N      = cs.N;
			calib.dir    = cs.transect.dir();
			calib.center = cs.transect.centre;
			% TODO unsafe for reg
			calib.dw     = cs.grid_n.dx1;
			if (~isempty(cs.grid_nz))
				calib.dz     = cs.grid_nz.dx2;
			end
			firsttime = false;
		end

		n = n+nt;
		clear cs
		end % if filename_C(idx)
	end % for idx

	% flux and depth averaged concentration
	%calib.qs_tn = calib.var_tn('sediment_discharge')./calib.q_tn;
	calib.ssc_tn = calib.sediment_discharge_tn./calib.q_tn;


%	calib.vbar  = vbar;
%	calib.Un    = squeeze(nanmean(calib.U,2));
%	calib.angle = 90+rad2deg(atan2(calib.ubar,calib.vbar));

	% average bed level
	%calib.lH = nanmean(calib.zb.A(:));

	% average over all campaigns
	calib.ln_z0.mean = nanmean(calib.ln_z0.A.val,2);
	calib.zb.median  = nanmedian(calib.zb.A,2);
	calib.zb.q16  = quantile(calib.zb.A',0.16)';
	calib.zb.q84  = quantile(calib.zb.A',0.84)';
	calib.zb.mean = mean(calib.zb.A,2);
	calib.zb.std  = std(calib.zb.A,[],2);
%	calib.zb.wmean  = sum(calib.zb.A.*calib.m1,2)./sum(calib.m1,2);


	% average bed level
	calib.zb0 = calib.zs0 - calib.radius;

	% test
	if (any(any(~isreal(calib.ln_z0.A.val))))
		warning('imaginary values in z0');
		calib.ln_z_0.A.val = real(calib.ln_z0.A.val);
	end

		calib.ln_z0.A.val = real(calib.ln_z0.A.val);
		calib.u_s.A.val = real(calib.u_s.A.val);
		calib.perturbation.A.val = real(calib.perturbation.A.val);
	
end % load_vadcp_discharge

