% Tue 18 Sep 16:38:09 CEST 2018
%% velocity near the bed
function [ubed,usur] = velocity_near_bed(obj,field)
	if (nargin()<2)
		field = 'earth';
	end
	S    = obj.bin.S;
	last = obj.ens.last;
	
	% allocate memory
	ubs = zeros(2,obj.ens.n,class(obj.velocity.(field)));	
	vbs = zeros(2,obj.ens.n,class(obj.velocity.(field)));	

	% for each ensemble
	for idx=1:obj.ens.n
		s = S(1:last(idx),idx);
		% bot, top
		A = [s,1-s];
		ubs(:,idx) = A \ obj.velocity.(field)(1:last(idx),idx,1);
		vbs(:,idx) = A \ obj.velocity.(field)(1:last(idx),idx,2);
	end % for idx
	obj.ens.velocity.bed.(field)     = [ubs(2,:).', vbs(2,:).'];
	obj.ens.velocity.surface.(field) = [ubs(1,:).', vbs(1,:).'];
end  % velocity_near_bed

