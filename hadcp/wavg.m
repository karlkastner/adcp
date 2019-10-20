% 2015-01-22 16:14:34.638913204 +0100
% Karl Kastner, Berlin
%% weighted average ?
function [m, w] = wavg(c,u,s2c,s2u)
	for idx=1:size(c,2)
		p(:,idx) = wopt(c(:,idx),u(:,idx),s2c,s2u);
	end
	w = p.*c;
	m = cvec(sum(w.*u));
end


