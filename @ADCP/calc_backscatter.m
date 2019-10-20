% Sun Jun  8 13:58:23 WIB 2014
% Thu Dec  4 18:04:11 CET 2014
% Karl Kastner, Berlin
%
%% backscatter from echo intensity
% c.f. Deines; Wall 2006
%
function obj = calc_backscatter(obj)

	% filter the backscatter
	echo = obj.echo;
	echo(0 == echo) = NaN;
%	echo = nanmean(echo,3);
	% number of filtering steps, 1 for no-reweighing
	nf = 1;
	for idx=1:size(echo,3)
		echo(:,:,idx) = filteriir(echo(:,:,idx),obj.ens.H,obj.ens.last,obj.pf_bs,nf);
	end
	% write back
	obj.echo = echo;

%	% prepare filter
%	if (~isempty(obj.filter.mode))
%		obj.filter.n = round(obj.filter.t/((obj.time(2)-obj.time(1))*(3600*24)));
%		% prepare the filter window
%		obj.filter.win = gausswin(obj.nf);
%		% normalise window
%		obj.filter.win = obj.filter.win/sum(obj.filter.win);
%	end
	   	
	% near field correction factor
	% TODO, this has to be different for each hadcp beam
	psi = obj.near_field_correction();

	% array lengths
	nr     = size(obj.R,1);
	nt     = length(obj.time);

	% repeat arrays to equal size
	alpha_w = obj.sound_absorption_water();
	aa      = repmat(alpha_w', [nr, 1]);

	% TODO not true for HADCP
	RR     = obj.bin.R;
	T_C    = obj.transducer_temperature_C();
	T_K    = Constant.celsius_to_kelvin(T_C);
	TT     = repmat(rvec(T_K), [nr, 1]);
	PP     = repmat(rvec(obj.power_W), [nr, 1]);
	PsiPsi = zeros([max(obj.nbins),obj.nt],'single');
	L      = zeros([max(obj.nbins),obj.nt],'single');
	lngthtranspulse = obj.lngthtranspulse;
	for idx=1:obj.lastFile
		fdx = (idx == obj.dat.FileNumber);
		PsiPsi(:,fdx,:) = repmat(psi(:,idx),1,sum(fdx));
		L(:,fdx,:) = lngthtranspulse(idx);
	end

	% noise removal
	E = obj.echo;
	switch (obj.noiseflag)
	case {'log'}
		for idx=1:obj.nbeams
			E(:,:,idx) = E(:,:,idx) - obj.En(idx);
		end
	case {'lin'}
		% gostiaux van haren 2010
		% does not make a difference as ong as E - En > 10
		for idx=1:obj.nbeams
			% quick fix
			emin = -Inf;
			En(:,:,idx) = 10*log10(max(emin, ...
			  10.^(0.1*obj.Kc(idx)*E(:,:,idx)) ...
			- 10.^(0.1*obj.Kc(idx)*obj.En(idx)) ));
		end
	otherwise
		error('here');
	end

	% backscatter in dB, uncorrected for sediment attenuation
	% TODO correct Kc for temperature: 
	% according to shields-2010, kc~127.3/T_K
	% kc = kc_20C*c2k(20C)/T_k = 293.16/Tk*kc_20C
	S_a = zeros(max(obj.nbins),obj.nt,'single');
	for idx=1:obj.nbeams
		% TODO, no magic numbers
		if (isa(obj,'HADCP') && 3 == idx)
			% TODO also for psi
			% disp('hadcp found');
			RR_ = cos(obj.beamangle_rad(1))*RR;
			L_  = cos(obj.beamangle_rad(1))*L;
		else
			RR_ = RR;
			L_ = L;
		end
		S_a(:,:,idx) = ...
		          obj.Kc(idx)*E(:,:,idx) ... % counts to db
			+ 2*aa.*RR ...	
			+ 10*log10(PsiPsi.^2.*TT.*RR.^2./(L.*PP)) ...
			+ obj.C;
	end

	obj.backscatter_log10 = S_a;

%	% filter the backscatter values in geometric space
%	if (strcmp(obj.filter.mode,'log'))
%		% TODO this is problematic when depth changes
%		obj.S_a = conv2_man(1,obj.filter.win,obj.S_a,'same');
%	end

	% backscatter in arithmetic space
	obj.backscatter  = 10.^(0.1*S_a);

	% filter the backscatter in arithmetic space
	if (strcmp(obj.filter.mode,'lin'))
		obj.beta = conv2_man(1,obj.filter.win,obj.beta,'same');
	end

	if (0)
	% integral of backscatter along the beam
	% TODO, account for differently sized first cell
	% TODO dr varies for the HADCP beams
	dr = obj.dr;
	obj.Ibeta = zeros(size(obj.backscatter));
	for idx=1:obj.lastFile
		fdx = (idx == obj.dat.FileNumber);
		obj.Ibeta(:,fdx,:) = dr(idx)*cumsum(obj.backscatter(:,fdx,:),1);
	end
	end

end % calc_backscatter

