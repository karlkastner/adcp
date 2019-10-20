% Sun Jan  5 00:06:42 WIB 2014
% Fri 26 Jan 16:46:59 CET 2018
% Karl Kastner, Berlin
%
%% process VADCP data
function obj = process(obj,transect)

	% TODO preallocate values in vel profile and backscatter as NaN (!)
	for cdx=1:length(transect)
		ensmask = logical(obj.ens.tid(:,cdx));

	        % rotate the velocity vectors to local flow direction
		obj.to_cs(transect(cdx).dir,ensmask);
		obj.depth_average_velocity('earth',ensmask); 
		obj.depth_average_velocity('cs',ensmask); 

		obj.to(obj.shear_stress_field);
		obj.fit_velocity_profile(obj.shear_stress_field,ensmask);

		obj.velocity_near_bed();

		obj.calc_backscatter();
		obj.calc_ssc();
		
		obj.depth_integrate_sediment_discharge('cs',ensmask);

		obj.fit_sediment_concentration_profile(ensmask);
	end % cdx
end % process

