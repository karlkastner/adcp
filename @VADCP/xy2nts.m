% Sat Jan 25 13:35:28 WIB 2014
% Karl Kastner, Berlin
%
%% project coordinates onto a single cross section and assign them nz-coordinates at a single cross section
%% TODO this should be part of transect
%
function [obj] = xy2nts(obj,cdx,centre,dir,dwidth,T_max)
	flag = false;

	% average coordinate of all beams
	[obj.N(:,cdx), obj.T(:,cdx)]  = xy2nt(obj.X, obj.Y, centre, dir, flag);

	% coordinate of each beam
	[obj.N4(:,:,cdx), obj.T4(:,:,cdx)] = xy2nt(obj.ens.X4,obj.ens.Y4, centre, dir, flag);

	convex  = abs(obj.N(:,cdx)) <= 0.5*dwidth;
	inrange = abs(obj.T(:,cdx)) < T_max;

	% test, if sample is within the transect in spanwise direction
	%       close enough in the streamwise direction
	%	no value has yet been assigned or this cross section is closer
	% TODO  the latest implementation allows samples to belong to multiple cross sections
	fdx_ =  convex & inrange;
	      %... % & ((0 == tid) ...
	      %... | abs(T_) < abs(T));

%	obj.ens.tid(1:obj.ens.n,cdx) = false;
	obj.ens.tid(fdx_,cdx) = true;
end % xy2nts

