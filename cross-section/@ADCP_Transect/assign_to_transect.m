% Sun Jan  5 00:06:42 WIB 2014
% Karl Kastner, Berlin
%
% Fri 26 Jan 16:46:59 CET 2018
%
%% assign ensemble to respective transects
%% this has a side-effect (writes to) the adcp object,
%% but values of induvidial cross sections remain unaffected by each other
%
function obj = assign_to_transect(obj)
        % determine the cross section
        % and rotate the velocity values with respect to cross section orientation
        obj.transect.fit();

	% filter velocity and backscatter data
	% TODO this is deprecated
	%f = 1;
	%f (nf > 1)
	%vadcp.velocity.earth = medfilt1(double(vadcp.velocity.earth),nf,[],2);
	%vadcp.beta      = medfilt1(double(vadcp.beta),nf,[],2);
	%nd

        % transform ensemble coordinates to cross section coordinates (n-t-s)
        % and determine to which cross section each ensemble belong
	obj.xy2nts();

        % identify into individual crossings (transects)
        adcp = obj.split_transect();    
end % Transect/preprocess_adcp()

