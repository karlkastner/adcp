% Thu 20 Sep 08:50:03 CEST 2018
%% transform velocity to given reference
function to(obj,field)
	switch (field)
	case {'abs'} % magnitude and direction
		obj.to_abs();
	case {'earth'}
		obj.to_earth();
	case {'beam'}
		obj.to_beam();
	case {'instrument'}
		obj.to_beam();
	case {'ship'}
		obj.to_ship();
	case {'cs'} % w/r to cross section
		obj.to_cs();
	case {'sw'} % stream and spanwise
		obj.to_sw();
	otherwise
		error('here');
	end
end
