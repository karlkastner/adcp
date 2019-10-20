% Thu Jul 17 21:26:14 WIB 2014
% Karl Kastner, Berlin
%
%% depth of bin, distance between water surface (z_s) and (z_i)
%%
%% Ds = z_s - z_bin
%%
%% does not correct for tilts
%
% TODO this has to be defined for every ensemble
function Ds = Ds(obj)
	if (~obj.is_facing_upward)
		Ds = obj.d_transducer + obj.Dt;
	else
		%Ds = bsxfun(@obj.Dt
	end
end

