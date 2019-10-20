% 2013-10-16 21:51:04 +0700
% Karl Kastner, Berlin

%% evaluate the objective function at the selected points

function [XX YY f g] = bsgrid(obj,n,m)
	if (nargin() < 2)
		n = 200;
	end
	if (nargin() < 3)
		m = 4;
	end
	x0 = obj.param.mean;
	X = zeros(1,n);
	Y = zeros(1,n);
%	for idx=1:length(x0)
		X(:) = linspace(0,x0(1)*m,n);
		Y(:) = linspace(0,x0(2)*m,n);
%	end
	nc = length(obj.M_ref);
	S = obj.param.cov;
	iS = inv(S);
	mu = obj.param.mean';
	for idx=1:n
	for jdx=1:n
		M_est = obj.func([X(idx),Y(jdx)], ...
			obj.I_f_1, obj.I_f_2, obj.beta_f_1, obj.beta_f_2, obj.PC_ref, obj.r_ref);
		e = M_est - obj.M_ref;
		f(idx,jdx) = sqrt(e'*e/nc);
		%g(idx,jdx) = obj.param.cov*[X(idx);
		d = (mu - [X(idx); Y(jdx)]);
		g(idx,jdx) = obj.rmse + sqrt(nc*d'*iS*d);
		XX(idx,jdx) = X(idx);
		YY(idx,jdx) = Y(jdx);
		%g(idx,jdx) = obj.rmse + sqrt(nc*d'*iS*d);
	end
	end
end

