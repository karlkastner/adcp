% Thu Mar  5 17:36:47 CET 2015
% Karl Kastner, Berlin
%
%% number of ensembles
%
function nens = nens(obj)
	% TODO should be bin.n and ens.n
	nens = size(obj.dat.VEL,2);
end

