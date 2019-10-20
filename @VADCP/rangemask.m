% Thu Jul 17 19:31:07 WIB 2014
% Karl Kastner, Berlin
%% mask all bins in range
function [rangemask] = rangemask(obj,varargin)
	if (isempty(obj.rangemask_))
		if (~isempty(obj.btrange))
		beamangle_rad = obj.beamangle_rad;
		btrange = obj.btrange;
		Dt      = obj.Dt;

		% mask samples that are unaffected by bottom side lobes
		% side lobes propagate on the shortest way to the bottom and back,
		% i.e. vertically dowm, therefore beams, which are slanted,
		% are cut off at a radius equal the minimum depth
		min_btrange = min(btrange,[],2);
		nFile = size(Dt,2);
		for idx=1:nFile
			fdx = (idx == obj.dat.FileNumber);
			rmask(:,fdx) = bsxfun(@lt,Dt(:,idx),cos(beamangle_rad(1))*min_btrange(fdx)');
		end
		obj.rangemask_ = rmask;
		else
			obj.rangemask_ = true(obj.bin.n,obj.ens.n);
		end
	end
	if (nargin() < 2 )
		rangemask = obj.rangemask_;
	else
		rangemask = obj.rangemask_(:,varargin{:});
	end
end

