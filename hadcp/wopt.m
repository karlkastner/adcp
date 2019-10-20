% 2015-01-22 16:14:34.638913204 +0100
% Karl Kastner, Berlin
%% optimal weights for averaging (lumped) velocities that are each associated
%% with error variance s2
function p = wopt(c,u,s2c,s2u)
	n = size(u,1);
	b = ones(n-1,1)*(c(end)^2*s2u(end) + u(end)^2*s2c(end));
	A = ones(n-1)*(c(end)^2*s2u(end) + u(end)^2*s2c(end)) ...
	    + diag( u(1:end-1).^2.*s2c(1:end-1) + c(1:end-1).^2.*s2u(1:end-1) );
	p        = A\b;
	p(end+1) = 1-sum(p);
end

