% Thu 11 Jul 11:38:07 CEST 2019
%% plot the velocity components
% TODO time vs ensemble
function plot_velocity_components(obj,field)
	if (nargin()<2)
		field = 'earth';
	end
	v=obj.velocity.(field);
	m = obj.mask;
	l={'u','v','w','e'};
	for idx=1:4;
%		figure(1);
		subplot(2,2,idx);
		cla();
		v_=v(:,:,idx);
		v_(~m)=NaN;
		v(:,:,idx)=v_;
		%imagesc(v_);
		surface(v_,'edgecolor','none');
		axis ij
		c = colorbar;
		xlabel(c,'m/s');
		title(l{idx});
		xlabel('ensemble');
		ylabel('bin');
	end % for 
end % function

