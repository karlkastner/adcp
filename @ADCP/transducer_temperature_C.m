% Mi 10. Sep 15:29:20 CEST 2014
% Mon Jul 28 16:39:42 WIB 2014
% Karl Kastner, Berlin
%
%% convert raw transducer temperature to SI units [Celsius]
%% T   : (1,nt)  water temperature
% Note: used to be Kelvin, but this was wrong
function t_C = transducer_temperature_C(obj)
		% raw temperature
		t_raw = single(obj.dat.ADC(:,6))*256;
		% filter invalid values
		t_raw(t_raw == 0) = NaN;
		% temperature offset
		t_offset_C = obj.temperature_offset_C();
		% TODO multifile
		t_offset_C = t_offset_C(1);

		% real-time temperature of the transducer (deg C)
		c   = obj.TEMP_ADC;

		% c.f. WorkHorse H-ADCP Operation Manual, May 2015
                t_C = t_offset_C + (((c(4).*t_raw + c(3)).*t_raw + c(2)).*t_raw + c(1));
	    	% temperature_K = (ADCP.TEMPOFFSET_K + ADCP.TEMPSLOPE_K*single(temperature_raw))';
end % transducer_temperature_C

