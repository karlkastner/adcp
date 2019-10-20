% Sat Nov  2 17:49:42 UTC 2013
% Karl Kästner, Berlin
%% new fiel correction of the acoustic backscatter
%% c.f. wall 2006
%% Psi : (nr,1)  near field correction factor
function [psi, obj] = near_field_correction(obj)
	if (obj.psiflag)
		% TODO, this is an approximation, but the sound velocity is not so much variable
		cw     = obj.sound_velocity();
		cw     = nanmean(cw);
		lambda = cw/obj.FREQ_HZ;
		% transceiver radius
		rn     = pi*obj.AT_M^2/lambda;
	
		% for each file
		R = obj.R;
		for idx=1:obj.lastFile
			z(:,idx) = R(:,idx)/rn;
		end
		psi =    ( 1 + 1.35*z + (2.5*z).^3.2 ) ...
                      ./ ( 1.35*z + (2.5*z).^3.2 );
	
	else
		psi = ones(max(obj.nbins),'single');
	end
% klammer falsch in der publication -> offensichtlich nicht
	% obj.psi = ( 1 + 1.35*z + 2.5*z.^3.2 ) ./ ( 1.35*z + 2.5*z.^3.2 );
	%obj.psi = obj.psi.^0.5;
	%obj.psi = ones(size(obj.psi));
	%obj.psi = sqrt(obj.psi);
end % method near_field_correction()

