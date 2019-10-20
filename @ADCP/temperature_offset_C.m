% Tue  8 May 08:47:58 CEST 2018
% Temp Sens Offset (PS0 results)
%% instrument specific temperature offset
function t_offset_C = temperature_offset_C(obj)
	serial = obj.serial;
	for idx=1:length(serial)
	% TODO, load from csv
	switch (serial(idx))
	case {2525}
		t_offset_C(idx) = -0.30;
	case {9534}
		t_offset_C(idx) = -0.23;
	case {2549}
		t_offset_C(idx) = -0.27;
	case {9964}
		t_offset_C(idx) = -0.09;
	case {9974}
		t_offset_C(idx) = -0.15;
	case {1094}
		t_offset_C(idx) = 0;	% not given by firmware for stream pro
	otherwise
		t_offset_C(idx) = 0;
		warning('no temperature offset for instrument found, using 0');
	end % switch
	end % for
end

