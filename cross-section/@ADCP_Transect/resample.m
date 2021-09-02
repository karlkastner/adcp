% Wed 21 Oct 16:10:23 +08 2020
% TODO make field a cell (multiple fields)
% TODO option for averaging or not
function s = resample(obj,field,vadcp,siz)
	% oversampling
	fn = 2;

	%NN = vadcp.bin.N;
	SS = vadcp.bin.S;
	val = vadcp.(field);
	val = mean(val,3);

	% invalidate last bin
	ind=sub2ind(size(val),max(vadcp.last_bin,1),(1:size(val,2))');
	val(ind)=NaN;

	msk = vadcp.rangemask;
	SS  = vadcp.bin.S;

	% for each cross-section
	for odx=1:length(obj)
		cdx = obj(odx).tdx;
		% number of cells along vertical
		if (nargin > 3)
			ns = siz(odx).n(1);
			nn = siz(odx).n(2);
		else
		% TODO make dynamic along N, depending on local depth
		ns  = fn*max(vadcp.ens.last(vadcp.ens.tid(:,odx)));
		% number of cells along horizontal
		nn  = fn*(median(obj(odx).last-obj(odx).first)+1);
		end
		s(odx).N     = innerspace(-0.5*obj(odx).dwidth,0.5*obj(odx).dwidth,nn);
		s(odx).bin.N = repmat(s(odx).N,ns,1),
		s(odx).bin.S = repmat(innerspace(1,0,ns)',1,nn);
		% for each transect
		nt = 0;

		fail = false(size(obj(odx).first));
		for idx=1:length(obj(odx).first)
			N    = vadcp.N(obj(odx).first(idx):obj(odx).last(idx),cdx);
			NN_  = repmat(rvec(N),size(val,1),1);
			SS_  = SS(:,obj(odx).first(idx):obj(odx).last(idx));
			val_ = val(:,obj(odx).first(idx):obj(odx).last(idx));
			msk_ = ( msk(:,obj(odx).first(idx):obj(odx).last(idx)) ...
				& isfinite(val_) & isfinite(NN_));
			msk_ = flat(msk_);
			if (sum(msk_(:))>2)
			nt = nt+1;
			sc = 7;
			interp       = TriScatteredInterp(double([NN_(msk_),sc*SS_(msk_)]),double(val_(msk_)));
				
			s(odx).(field)(:,:,nt) = reshape(interp([flat(s(odx).bin.N),sc*flat(s(odx).bin.S)]),ns,nn);
			
			% depth
			H       = vadcp.ens.H(obj(odx).first(idx):obj(odx).last(idx));
			[N,sdx] = sort(N);
			H = H(sdx);
			ddx = [true; N(2:end) ~= N(1:end-1)] & isfinite(H) & isfinite(N);
			s(odx).H(:,nt) = interp1(N(ddx),H(ddx),s(odx).N,'linear');

			else
				fail(idx) = true;
			end
		end % for idx (each crossing)
		s(odx).field(:,:,fail) = NaN;
		s(odx).H(:,fail) = NaN;
	        s(odx).bin.Z = -s(odx).bin.S.*nanmedian(s(odx).H,2)';
	end % for odx (each cross-section)

if (0)
		% resample along vertical
		% TODO should be function in vadcp
		fdx = find(rvec(v.ens.tid(:,odx)));
		for jdx=1:length(fdx)
			zd = 1:v.ens.last(fdx);
			bs_(:,jdx) = bs(zd,fdx(jdx));
			% TODO use distance to transducer (!)
			bs__(:,jdx) = interp1(zd(flag)/(v.ens.last(fdx)+1), ...
					      bs(zd(flag),fdx(jdx)), ...
					      (1:ns)/(ns+1), 'linear');
		end
		% TODO should be a function of ad
		% resample along the horizontal
		bs_s(odx).N = linspace(-obj(idx).dwidth/2,+obj(idx).dwidth);
		for jdx=1:length(obj(odx).first)
			N       = v.N(obj(odx).first(jdx):obj(idx).last(jdx));
			[N,sdx] = sort(N);
			bs_s(odx).echo(:,:,jdx) = interp1(N,bs__(:,sdx),bs_s(idx).N,'linear');
			% TODO, depth
			H       = v.ens.H(obj(odx).first(jdx):obj(idx).last(jdx));
			H = H(sdx);
			N = cvec(N);
			ddx = [true; N(2:end) ~= N(1:end-1)] & isfinite(H);
			%N = cvec(N)+(0:length(N))'*sqrt(eps)*max(abs(N));
			bs_s(odx).H(:,jdx) = interp1(N(ddx),H(ddx),'linear');
		end 
end
end % resample

