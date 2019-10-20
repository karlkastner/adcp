% So 30. Aug 09:05:03 ECT 2015
% Karl Kastner, Berlin
%% gui user selection of cross-section end points
function [x0, y0] = define_transect(adcp,m)

	if (nargin() < 2)
		m=1;
	end
	
	[X Y]=utmADCP(adcp);
	
	figure(1)
	clf();
	scatter(X,Y,1);
	axis equal;
	hold on;
	
	x0=[];
	y0 =[];
	idx = 0;
	while(true)
		idx = idx+1;
		for jdx=1:2
			printf('Transect %d end point %d\n',idx,jdx);
			[x0(idx,jdx) y0(idx,jdx)] = get_coordinates(1);
			if (strcmp('alt',get(gcf,'Selectiontype')))
				x0 = x0(1:idx-1,:);
				y0 = y0(1:idx-1,:);	
				return;
			end
			hold on
			plot(x0(idx,jdx),y0(idx,jdx),'ro');
		end
		plot(x0(idx,:),y0(idx,:),'r-','linewidth',2);
	end
	%save('../dat/amazon-endpoints.mat');
end

