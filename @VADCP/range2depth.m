% Mon  4 Jun 10:45:53 CEST 2018
%% depth below transducer for individual bins of the beams
function D = range2depth(obj,R,roll,pitch,beamangle_rad,invert)
	if (nargin() < 6)
		invert = false;
	end 
	if (isempty(R))
		D = [];
		return;
	end
	if (1)
		% c.f. A. Woodgate and A.E. Holroyd, 2011
		% formula had to be adapted for downward looking ADCPS (e.q. of ADCP Coordinate Transformation, p. 19)
		% adapt roll for downward looking deployment
	%	roll = -(roll-deg2rad(180));
		% adapt pitch as pitch and roll not measured by gimbals (effecting each other)
	%	pitch = atan(tan(pitch).*cos(roll));
		cos_a = sin(beamangle_rad)*[-sin(roll(:)), +sin(roll(:)), +sin(pitch(:)), -sin(pitch(:))] ...
			+ repmat(cos(beamangle_rad)*sqrt(1 - sin(roll(:)).^2 - sin(pitch(:)).^2), 1, 4);
		if (~invert)
			if(size(R,2) == size(cos_a,2))
				D = R.*(cos_a./cos(beamangle_rad));
			else
				for idx=1:size(cos_a,2)
					D(:,:,idx) = R*(cos_a(:,idx)./cos(beamangle_rad)).';
				end
			end
		else
			D = bsxfun(@times, cos(beamangle_rad)./cos_a, R);
		end
	else
	% TODO, for some reason rotation does _not_ work	
		for idx=1:size(btrange,1)
		RP = [	1 0 0
	      		0 cos(pitch(idx)) -sin(pitch(idx))
			0 sin(pitch(idx))  cos(pitch(idx))];
		RR = [ cos(roll(idx))  0  sin(roll(idx))
		       0       1       0
		       -sin(roll(idx)) 0  cos(roll(idx)) ];
		v = [	          0,          0, tan(beamangle_rad), -tan(beamangle_rad);
			-tan(beamangle_rad), tan(beamangle_rad),          0,           0;
			          1,          1,          1,           1 ];
		v = RP*RR*v*diag(btrange(idx,:));
		%	v = RR*RP*v;
		btrange(idx,:) = v(3,:);
	end
end

