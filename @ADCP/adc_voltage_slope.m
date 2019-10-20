% Tue  8 May 09:17:37 CEST 2018
%% instrument series specific conversion factors for voltage
%% c.f. WorkHorse Commands and Output Data Format, March 2014
%% c.f. XMT Voltage and Current Channels
%% originally undoccumented by RDI, and taken from Shields 2010
function slope_V = adc_voltage_slope(obj)
	switch (round(obj.FREQ_HZ/1e3))
	case {77}
		slope_V = 2092719/1000000;
	case {154,150}
		slope_V = 592157/1000000;
	case {307,300}
		slope_V = 592157/1000000;
	case {614,600}
		slope_V = 380667/1000000;
	case {1229,1200}
		slope_V = 253765/1000000;
	case {2458,2400}
		slope_V = 253765/1000000;
	otherwise
		slope_V = NaN;
		warning(sprintf('could not  load adc voltage slope for frequency %f',obj.FREQ_HZ));
	end
end

