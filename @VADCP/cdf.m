% Thu 11 Jul 11:05:10 CEST 2019
%% compute and plot cumulative distribution (cdf) of the velocity components
% TODO, make resampling and plotting optional
function [x,p] = cdf(obj,field)
	if (nargin()<2)
		field = 'earth';
	end
	v = obj.velocity.(field);
	m = obj.mask;
	m = repmat(m,[1,1,4]);
	v(~m) = NaN;
	s = size(v);
	n = prod(s(1:2));
 	v=reshape(v,[n,4]);
	v = sort(v);
	nf = sum(isfinite(v));
	x = (1:n)'*(1./nf);
	% resample
	ni = 100;
	xi = innerspace(0,1,ni)';
	vi = zeros(ni,4);
	for idx=1:4
		vi(:,idx) = interp1(x(1:nf(idx)),v(1:nf(idx),idx),xi);
	end
	plot(vi,xi);
	if (0) % pdf
		vi_ = csmooth(bsxfun(@times,diff(xi),1./diff(vi)),10);
		plot(mid(vi),vi_);
	end
	l={'u','v','w','e'};
	legend('location','northwest',l)
	xlabel('m/s')
	ylabel('cdf');
	grid on
end

