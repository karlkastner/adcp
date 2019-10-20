% Wed Oct 29 10:43:55 CET 2014
% Karl Kastner, Berlin
%
%% process backscatter, i.e. fit to cross-section grid from bin-values
%
% TODO, this function is superfluous
% vgrid is transposed with respect to the order of dimensions in the ADCP structure,
% so that the first dimension coincides with the bottom grid
function obj = process_backscatter(obj, adcp, ensmask)
	if (nargin()<3)
		ensmask = logical(adcp.ens.tid(:,obj.cdx));
	end
%	switch (lower(obj.vdimension))
%	case {'n','tn'}
		obj.process_backscatter_tn(adcp,ensmask);
%	case {'nz','tnz'}
%		obj.process_backscatter_tn(adcp,ensmask);
		obj.process_backscatter_tnz(adcp,ensmask);
%	otherwise
%		error('invalid vmethod');
%	end % switch
end % function CrossSection::process_backscatter()

