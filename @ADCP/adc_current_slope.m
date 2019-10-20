% Tue  8 May 09:16:58 CEST 2018
%% instrument type specific slope for converting raw current to Ampere
%% c.f WorkHorse Commands and Output Data Format, March 2014
%% c.f. XMT Voltage and Current Channels
%% originally undoccumented by RDI, and taken from Shields 2010
function slope_A = adc_current_slope(obj)
	switch (round(obj.FREQ_HZ/1e3))
	case {77}
		slope_A = 43838/1000000;
	case {154,150}
		slope_A = 11451/1000000;
	case {307,300}
		slope_A = 11451/1000000;
	case {614,600}
		slope_A = 11451/1000000;
	case {1229,1200}
		slope_A = 11451/1000000;
	case {2458,2400}
		slope_A = 11451/1000000;
	otherwise
		slope_A = NaN;
		warning('could not  load adc current slope');
	end
end

