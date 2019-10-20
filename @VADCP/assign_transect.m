% Sun Jan  5 00:06:42 WIB 2014
% Karl Kastner, Berlin
% Fri 26 Jan 16:46:59 CET 2018
%
%% assign transect index to ensembles
function obj = assign_transect(obj,transect)
	nt = length(transect);
	obj.ens.tid = false(obj.ens.n,nt);
% TODO N,T,N4,T4 should be part of ens
	obj.N   = zeros(obj.ens.n,nt);
	obj.T   = zeros(obj.ens.n,nt);
	obj.N4  = zeros(obj.ens.n,4,nt);
	obj.T4  = zeros(obj.ens.n,4,nt);

        % transform ensemble coordinates to cross section coordinates (n-t-s)
        % and determine to which cross section each ensemble belong
	for cdx = 1:nt
		obj.xy2nts(cdx,transect(cdx).centre,transect(cdx).dir,transect(cdx).dwidth,transect(cdx).T_max);
	end
end % assign_transects

