% Sat 14 Jan 23:13:20 CET 2017
%% discharge division ratio
function d = discharge_division(Q)
	d = (Q(:,1)-Q(:,2))./(Q(:,1)+Q(:,2));
end

