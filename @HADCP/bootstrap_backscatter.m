% Sun Jun  8 14:09:03 WIB 2014
% Karl Kastner, Berlin
%% bootstrap uncertainty of the backscatter parameters
function obj = bootstrap_backscatter(obj, n)
	t_ref = obj.t_ref;
	M_ref = obj.M_ref;
	M_reg = obj.M_reg;
	param = obj.param;
	m = length(M_ref);
	for idx=1:n
		% obj.t_ref = t_ref(p(:,idx));
		% obj.M_ref = M_ref(p(:,idx));
		% do generate permutationint the loop, as otherwise out of memory may occur
		p = randi(n,m,1);
		obj = feval(obj.regC{obj.model},obj,obj.carg{:}, p);
		obj.boot.param.B(idx,:) = obj.param.mean(:)';
		% todo, evaluation shoult take place not at random points
		%M = M + obj.M_reg;
	end
	obj.boot.param.std  = std(obj.boot.param.B);
	obj.boot.param.mean = mean(obj.boot.param.B);
	%obj.boot.M.mean     = M/n; 
	% restore
	obj.t_ref = t_ref;
	obj.M_ref = M_ref;
	obj.M_reg = M_reg;
	obj.param = param;
end % bootstrap backscatter

