% Sun Jul 13 12:16:49 WIB 2014
% Karl Kastner, Berlin
%% convert bytes of serial number into single number
function serial = convert_raw_serial_STATIC(serial)
		if (4 == size(serial,2))
			% construct the serial number from bytes
			% ADCP stores bytes in little endian order
			serial = cast_byte_to_integer(serial);
			%% big endian system
			%obj.serial =        256*uint32(serial(1,1)) ...
                        %             +        1*uint32(serial(1,2)) ...
                        %             + 16777216*uint32(serial(1,3)) ...
                        %             +    65536*uint32(serial(1,4));
		else
			warning('raw serial number does not expected format, skipping conversion');
		end
end % convert_raw_serial

