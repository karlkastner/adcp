% Fri Mar 27 21:11:08 CET 2015
%% homogenize the hadcp profile
function [bin, ens] = hadcp_homogenize_profile(time,U,rid,idmax,Ti,T,order,amporder,phaseorder)
	nbin = size(U,1);
	% maximum number of 
	dmax       = 4;
	robustflag = 1;

	fprintf('Short time fourier transform for each bin\n');
	phase = [];
	amp   = [];
	% for individual bins
	bin = struct();
	% ensemble averaged
	ens = struct();
	%U_  = NaN(size(U));
	% stft  = STFT();
	for idx=1:nbin
	%	disp(idx);
		stft(idx) = STFT('Ti',Ti,'T',T,'t0',time(1),'tend',time(end),'order',order);
		stft(idx).transform(time,U(idx,:).');

		[bin.U0(idx,:) void bin.Uf0(idx,:,:)] = stft(idx).itransform();

	%	U(idx,:) = stft(idx).itransform();
		phase(:,:,idx) = stft(idx).phase;
		amp(:,:,idx)   = stft(idx).amplitude;
	end

	fprintf('Removing outliers\n');
	% remove bad days, this usually coincides with data gaps and ends
	%isbad = [63 69:80 261 272 348:355 363:371 398 407:413 425:426 474:478 480];
	amax   = squeeze(max(max(amp(:,:,1:idmax),[],1),[],3));
	d      = cvec(abs(amax - nanmedian(amax)));
	isgood = d<4*nanmedian(d);
	
	phase(:,~isgood,:) = NaN;
	amp(:,~isgood,:)   = NaN;
	
	%for idx=1:nbin
	%	stft(idx).set_amplitude(amp(:,:,idx),phase(:,:,idx));
		%stft(idx).set_amplitude(amp(:,:,idx));
		%stft(idx).set_phase(phase(:,:,idx));
		% U_uc(idx,:) = stft(idx).itransform();
	%end
	fprintf('Removed %d intervals\n',sum(~isgood));

	fprintf('Phase lag profile');

	% phase difference
	dp.val    = bsxfun(@minus,phase,phase(:,:,rid));
	dp.val    = wrapToPi(dp.val);

	% linear order correction terms
	a1    = amp(1,:,rid);
	a2    = amp(2,:,rid);
	%a = amp(1,:,rid)+amp(2,:,rid);
	a1    = (a1 - nanmedian(a1))./diff(quantile(a1,[0.7 0.3]));
	a2    = (a2 - nanmedian(a2))./diff(quantile(a2,[0.7 0.3]));
	qa1   = quantile(a1,[0.1,0.9]);
	qa2   = quantile(a2,[0.1,0.9]);

	a = atan(amp(2,:,rid)./amp(1,:,rid));
	a = a-nanmedian(a);

	x = (1:size(amp,3))';

	%
	% average phase lag profile along profiling range
	% 
	nf      = size(amp,1);
	ns      = size(amp,2);
	nbin    = size(amp,3);
	dp.mean = [];
	% for each bin
	for idx=1:nbin
	   % for each frequency component
	   for jdx=2:nf
	        fdx = isfinite(cvec(a1)) ...
	                & isfinite(cvec(a)) ...
	                & isfinite(cvec(a2)) ...
	                & isfinite(cvec(dp.val(jdx,:,idx))) ...
	                & isfinite(cvec(dp.val(jdx,:,rid)));

		switch (phaseorder)
		case {0} % constant phase shift
			Xp = ones(ns,1);
		        Y = cvec(dp.val(jdx,fdx,idx));
%		        c = medianangle2(Y);
		        c = hodges_lehman(Y);

		        % TODO this is a generalised least squares problem
		        %X = ones(size(cvec(dp(jdx,fdx,idx))));
		        %c = mean(Y./X);
		        %c = X \ Y;
		        dp.mean(jdx,idx) = c(1);
		case {1}
			Xp  = [ones(ns,1) cvec(a)];
		        Y   = cvec(dp.val(jdx,fdx,idx));
			reg = Theil();
			c   = reg.fit(Xp(fdx,2),Y);
		        dp.mean(jdx,idx,:) = c;
		case {2}

		        % from local
			% TODO, this has to be done for sine and cosine !!!
			% but dp is small, so no large difference
		        %X = [ones(sum(fdx),1) cvec(a1(fdx)) cvec(a2(fdx))];
		        X = [ones(sum(fdx),1) cvec(a1(fdx)) cvec(a2(fdx))];

		        %X = [cvec(amp(jdx,fdx,idx)) cvec(a1(fdx).*amp(jdx,fdx,idx))];
		        Y = [cvec(amp(jdx,fdx,idx))];
		      	%X = ones(size(Y));
			% to reference
		        c = robustfit(X,Y,[],[],'off');
			%c = X \ Y;

		        % predict the mean profile
		        dp_mean(jdx,idx,:) = c;

		        % low flow and high flow profiles profiles
		        dp_mean_low(jdx,idx)  = [1 qa1(1)*1 qa2(1)*1]*c;
		        dp_mean_high(jdx,idx) = [1 qa1(2)*1 qa2(2)*1]*c;
		end % switch
	   end % for jdx
	end % for idx

	%
	% fit a function to the average phase lags to correct phase with 
	% respect to the channel centre and not just to reference bin
	% 
	fdx      = (x<=idmax);
	dp.mean_fitted = [];
	dp.coeff = [];
	for idx=2:nf
		%c = [0 0 0];
		% c2 cannot be taken into the exponent, due to sign
		f = @(c,x) c(1) + c(2)*exp(c(3)*x);
		%f = @(c,x) (c(1) + exp(c(2)*(x-c(3))));
		%c(1) = dp_mean(end,idx);
		c = [dp.mean(idx,idmax) -1 0];
		c = lsqnonlin(@(c) finite((f(c,x(fdx)) - dp.mean(idx,fdx).')),c);
		dp.coeff(idx,:)  = c;
		dp.mean_fitted(idx,:) = f(c,x);

		% TODO, this has actually to be done for each moment in time
%		c                    = lsqnonlin(@(c) w.*(f(c,x) - dp_mean1_low(idx,:)'),c);
%		C_low(idx,:)         = c;
%		dp_mean1_low_(idx,:) = f(c,x);
%
%		c                    = lsqnonlin(@(c) w.*(f(c,x) - dp_mean1_high(idx,:)'),c);
%		C_low(idx,:)         = c;
%		dp_mean1_high_(idx,:) = f(c,x);
	end

	% determine phase lag with respect to channel centre
	dp.mean_shifted        = bsxfun(@minus,dp.mean,dp.coeff(:,1));
	dp.mean_fitted_shifted = bsxfun(@minus,dp.mean_fitted,dp.coeff(:,1));

	% amplitude scales
	fprintf('Amplitude scale profile\n');
	
	% mean profiles, linear in mean flow velocity
	ra.val = bsxfun(@times,amp(:,:,rid),1./amp);
	ra.mean = [];
	for idx=1:nbin
	   for jdx=1:nf
		fdx = isfinite(cvec(a1)) ...
		        & isfinite(cvec(a2)) ...
			& isfinite(cvec(amp(jdx,:,idx))) ...
			& isfinite(cvec(amp(jdx,:,rid)));

		switch (amporder)
		case {0}
			% constant profile
			% TODO this is a generalised least squares problem
			if (robustflag)
				Xa = ones(ns,1); %cvec(amp(jdx,:,idx));
				Y  = cvec(ra.val(jdx,fdx,idx));
	
				% c = X \ Y;
				%c = median(Y./Xa(fdx));
				c = hodges_lehman(Y./Xa(fdx));
				ra.mean(jdx,idx) = c(1);
			else
				X = cvec(amp(jdx,fdx,idx));
				Xa = ones(ns,1);
				Y = cvec(amp(jdx,fdx,rid));
				ra.mean(jdx,idx) = X \ Y;
			end
		case {1}
			Xa  = [ones(ns,1) cvec(a)];
		        Y   = cvec(ra.val(jdx,fdx,idx));
			reg = Theil();
			c   = reg.fit(Xa(fdx,2),Y);
		        ra.mean(jdx,idx,:) = c;
		case {2}
			% from local
			X = [cvec(amp(jdx,fdx,idx)) cvec(a1(fdx).*amp(jdx,fdx,idx)) cvec(a2(fdx).*amp(jdx,fdx,idx))];
	
			%$X = [cvec(amp(jdx,fdx,idx)) cvec(a(fdx).*amp(jdx,fdx,idx))];
			% to reference
			%Y = [cvec(amp(jdx,fdx,rid))];
			c = X \ Y;
			%c = robustfit(X,Y,[],[],'off');
		%	c = robustfit(X,Y,'welsch',[],'off');
		%	c = quantreg(X,Y,0.5);
			c = theil2(X,Y,100*size(X,1));
	
			% predict the mean profile
			ra_mean1(jdx,idx,:) = cvec(c);
	
			% low flow and high flow profiles profiles
			ra_mean1_low(jdx,idx)  = [1 qa1(1)*1 qa2(1)*1]*c;
			ra_mean1_high(jdx,idx) = [1 qa1(2)*1 qa2(2)*1]*c;
		end % switch
	   end % for jdx
	end % for idx

	% fit a function through the average amplitude scales
	fdx      = (x<=idmax);
	ra.mean_fitted = [];
	ra.coeff = [];
	%C_low = [];
	%C_high = [];
	% correction for transverse gradient missed by HADCP
	cfac = 0.9;
	w = 450;
	for idx=1:nf
		c = [1 -1/w 0];
		%f = @(c,x) c(1)*(1-c(2)*exp(c(3)*x));
		f = @(c,x) 1./(c(1)*(1-exp(c(2)*(x-c(3)))));
		%c = [1 -1 0 0];
		%f = @(c,x) c(1)*(1-c(2)*exp(c(3)*x + c(4)*x.^2));
		c = lsqnonlin(@(c) (f(c,x(fdx)) - ra.mean(idx,fdx,1)'),c);
		%c = lsqnonlin(@(c) finite(bsxfun(@minus,f(c,x(fdx)),squeeze(ra(idx,:,fdx)).')),c);

		ra.coeff(idx,:)       = c;
		ra.mean_fitted(idx,:) = f(c,x);

		if (0)
		% high flow
		c                    = lsqnonlin(@(c) w.*(f(c,x) - ra_mean1_low(idx,:)'),c);
		C_low(idx,:)         = c;
		ra_mean1_low_(idx,:) = f(c,x);

		% low flow
		c            = lsqnonlin(@(c) w.*(f(c,x) - ra_mean1_high(idx,:)'),c);
		C_high(idx,:) = c;
		ra_mean1_high_(idx,:) = f(c,x);
		end
	end
	
	% 1/c(1) : scale with respect to the reference bin
	% c(2)   : decay rate of bank boundary layer
	% c(3)   : offset from bank, where velocity is zero
	c        = ra.coeff;

	% scale reference from bin rid to channel centre of infinite wide cross section
	ra.scale_inf                = c(:,1); %1./c(:,1);
	ra.mean_shifted_inf         = bsxfun(@times,ra.mean,ra.scale_inf);
	ra.mean_fitted_shifted_inf  = bsxfun(@times,ra.mean_fitted,ra.scale_inf);

	% scale reference from bin rid to cross sectional average
	% scale for finite width, bank offset ignored
	ra.scale = 1/cfac*c(:,1).*(1 - 2./(c(:,2).*w).*(exp(0.5*c(:,2)*w)-1));
	%1./c(:,1).*(1 - 2./(c(:,2).*w).*(1-exp(0.5*c(:,2)*w)));
	ra.mean_shifted  = bsxfun(@times,ra.mean,ra.scale);
	ra.mean_fitted_shifted  = bsxfun(@times,ra.mean_fitted,ra.scale);

	% shift with respect to infinity
	% U/ubar(0) = 2/w int)0^w/2 u_0 (1-exp(-n/L)) dn
	%scale_ = 1 - 2*c./w.*(1-exp(-0.5*c*w))

	% correct phase lag and amplitude scale
	siz=size(dp.mean);

	% phase and amplitude correction for each time step individually,
	% as profile is not constant in time for order>0
	ns = size(phase,2);
	for idx=1:nbin
	 for jdx=1:nf
		% predict instantaneous profile
		% TODO, shifted to cs average
		dp_(jdx,:,idx)  = Xp*squeeze(dp.mean(jdx,idx,:))';
		%if (~robustflag)
		%	Xa = cvec(amp(jdx,:,idx));
		%end
		ra_(jdx,:,idx) = Xa*squeeze(ra.mean_shifted(jdx,idx,:))';
         end
	end

	% fix for reference
	dp_(:,:,rid) = 0;
	ra_(:,:,rid) = 1; % TODO not true if shift is applied

	amp_c   = ra_.*amp;
	phase_c = phase - dp_;

	% arithmetic mean
	% cos and sin are linear combinations of the velocity, and have to be averaged, not phase and amplitude
	c = amp_c.*cos(phase_c);
	s = amp_c.*sin(phase_c);
%amp = nanmean(hypot(c(:,:,1:idmax),s(:,:,1:idmax)),3);
%u1 = nanmean(c(:,:,1:idmax)+s(:,:,1:idmax),3);
	c = nanmean(c(:,:,1:idmax),3).';
	s = nanmean(s(:,:,1:idmax),3).';
	avg.arith.amp      = hypot(c,s);
	avg.arith.phase    = atan2(s,c);
%u2 = avg.arith.amp.*(cos(avg.arith.phase) + sin(avg.arith.phase));
%narms(u1-u2,1)
%nanrms(amp'-avg.arith.amp)
%pause

	% harmonic mean estimate
	c = cos(phase_c).*amp;
	s = sin(phase_c).*amp;

	% combined scale + phase shift for real (cos part)
	%rc = ra_.*(cos(dp_)-sin(dp_));
	%rs = ra_.*(cos(dp_)+sin(dp_));
	valid = isfinite(s) & isfinite(c);
	%valid = isfinite(amp) & isfinite(phase);
	n     = sum(valid(:,:,1:idmax),3);
%	a     = amp;
%	a(~valid) = 0;
%	a     = sum(a(:,:,1:idmax),3)./n;
	c(~valid) = 0;
	s(~valid) = 0;
	c     = sum(c(:,:,1:idmax),3)./n;
	s     = sum(s(:,:,1:idmax),3)./n;
	ri    = 1./ra_;
	ri(~valid) = 0;
	r     = n./sum(ri(:,:,1:idmax),3);
	c     = (r.*c).';
	s     = (r.*s).';
	avg.harm.amp      = hypot(c,s);
	avg.harm.phase    = atan2(s,c);

	fprintf('Retransforming frequency components into velocities\n');

	% corrected velocities averaged along profile
	field_C={'arith','harm'};
	for idx=1:length(field_C)
		field = field_C{idx};
		avg.(field).stft = STFT('Ti',Ti,'T',T,'t0',time(1),'tend',time(end),'order',order);
		avg.(field).stft.set_amplitude(avg.(field).amp.',avg.(field).phase.');
		[avg.(field).U void avg.(field).U_] = avg.(field).stft.itransform(time);
	end

	% corrected velocity for each bin
	%bin = struct();
	bin.U = NaN(size(U));	
	for idx=1:nbin
		stft(idx).set_amplitude(amp_c(:,:,idx),phase_c(:,:,idx));
		%stft(idx).set_phase(c.phase(:,:,idx));
		%stft(idx).set_amplitude(c.amp(:,:,idx));
		[bin.U(idx,:) void bin.Uf(idx,:,:)] = stft(idx).itransform();
	end % for idx
	
	bin.time     = time;
	bin.amp      = amp;
	bin.phase    = phase;
	bin.amp_c    = amp_c;
	bin.phase_c  = phase_c;
	bin.ra       = ra;
	bin.dp       = dp;
	ens          = avg;
	ens.d        = d;
	ens.isgood   = isgood;
end % hadcp_homogenize_profile

