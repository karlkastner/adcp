% 2015-01-24 20:29:30.856568855 +0800
%% weighted average
function [m, w] = wavg(c,u,s2c,s2u)
%	for idx=1:size(c,2)
%		p(:,idx) = wopt(c(:,idx),u(:,idx),s2c,s2u);
%	end
%	w = p.*c;
%	m = cvec(sum(w.*u));
	m = mean(c.*u);
end

