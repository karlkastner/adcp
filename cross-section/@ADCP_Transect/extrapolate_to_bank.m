% 2017-01-19 19:04:07.529531875 +0100
%
%% extrapolate values to bank
% TODO do not take just first or last sample, but average over some range
function obj = extrapolate_to_bank(obj,adcp)
	% TODO, fit alpha for depending on the cross sections
	% TODO no magic number
	a = 0.5;

	N = adcp.ens.N;
	N = N(:,obj.tdx);

	% fetch
	first     = cvec(obj.first_valid);
	last      = cvec(obj.last_valid);
	sig       = cvec(obj.sig);

	area = obj.area;
	discharge = obj.discharge;


	% differentiate between left and right bank
	s        = N(first) < N(last);
	left     = first.*s +  last.*(~s);
	right    = last.*s  + first.*(~s);
	Nl       = N(left);
	Nr       = N(right);
	Hl       = adcp.ens.H(left);
	Hr       = adcp.ens.H(right);
	dir      = obj.dir;

	sig = s + (s-1);

	% projected discharge
	UV        = adcp.ens.velocity.earth;
	ql        = Hl.*(   dir(2).*UV(left,1) ...                  
                          - dir(1).*UV(left,2) ); 
	qr        = Hr.*(   dir(2).*UV(right,1) ...                  
                          - dir(1).*UV(right,2) ); 

	% discharge on left bank, exponential boundary layer
	dw = 0.5*obj.dwidth + Nl;
	Ql = sig.*ql.*(dw - (1-exp(-dw./Hl*a))/a);

	% linear extrapolation of area, 0 at boundary
	% TODO extrapolation with 30 deg bank
	area.left       = 0.5*Hl.*dw;

	% discharge on right bank
	dw = 0.5*obj.dwidth - Nr;
	Qr = sig.*qr.*(dw - (1-exp(-dw./Hr*a))/a);

	area.right      = 0.5*Hr.*dw;

	discharge.left = Ql;
	discharge.right = Qr;

	% totals
	area.total      = area.left + area.right + area.centre;
	discharge.total = discharge.left + discharge.right + discharge.centre;

	obj.discharge = discharge;
	obj.area      = area;
end % extrapolation_to_bank

