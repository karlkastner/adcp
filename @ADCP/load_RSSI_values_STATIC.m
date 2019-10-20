% Mi 10. Sep 15:37:39 CEST 2014
% Karl Kastner, Berlin
%% load instrument specific backscatter conversion parameters
function [Kc, En] = load_RSSI_values_STATIC(serial)
	Kc = NaN(1,4);
	En = NaN(1,4);

	w = what('@ADCP');                                                        
        p = w.path;    

	% TODO use classpath 
	t = readtable([p,filesep(),'rssi_values_calibrated.csv'],'delim',';');
	a = table2array(t);

	% add new adcps as rows
	% moved to rssi_csv.m
%	library.serial = uint32([2525;
%                                 2549
%                                 9964]);
%	library.Kc     = [ 0.3729  0.3739  0.3653  0.3693;
%                           0.3672  0.3721  0.3873  0.3723;
%                           0.4     0.4     0.4     0.4];
%	library.En     = [ 40         40         40         40;
%			   40         40         40         40;
%			   43.3351    43.7866    41.9715    42.7069];
	fdx = t.serial==serial(1);
	switch (sum(fdx))
	case {0}
		warning('no Kc value found, loading');
		Kc = 0.4*[1,1,1,1];
	case {1}
		Kc = a(fdx,2:5);
		En = a(fdx,6:9);
	otherwise
		error('dublicate entry');
	end
%	En = [40,40,40,40];
%	Kc = 0.4*[1,1,1,1];
%	for idx=1:size(library.serial,1)
%		%if (0 == sum(abs(library.serial(idx,:) - obj.serial)))
%		if (library.serial(idx) == serial)
%			disp('loaded Kc from library');
%			%Kc = library.Kc(idx,:);
%			En = library.En(idx,:);
%			break;
%		end % if
%	end % for idx
end % load_RSSI_values_STATIC

