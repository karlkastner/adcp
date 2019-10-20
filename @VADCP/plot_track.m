% Mon 28 May 08:20:36 CEST 2018
%
%% plot the boat track
function plot_track(adcp,transect,individual)
	if (nargin() < 2)
		transect = [];
	end
	if (nargin < 3)
		individual = false;
	end

	x0 = 0;
	y0 = 0;
	col = colormap();
        plot(adcp.X-x0,adcp.Y-y0,'.-');
        hold on;    
        for idx=1:length(transect)
                 %msk = logical(adcp.ens.tid(:,idx));                     
                 %plot(adcp.X(msk)-x0,adcp.Y(msk)-y0,'.','color',col(idx+1,:));
		first = transect(idx).first;
		last = transect(idx).last-1;
                 %msk = logical(adcp.ens.tid(:,idx));
		for kdx=1:length(first)
			if (individual)
				figure(kdx);
				clf();
        			plot(adcp.X-x0,adcp.Y-y0,'k.-');
			        hold on;    
			end
			msk=first(kdx):last(kdx);
			if (transect.isvalid(kdx))
        	        	plot(adcp.X(msk)-x0,adcp.Y(msk)-y0,'g.-'); %,'color',col(kdx+1,:));
			else
        	        	plot(adcp.X(msk)-x0,adcp.Y(msk)-y0,'r.-'); %,'color',col(kdx+1,:));
			end
			title(sprintf('ensemble %d - %d',first(kdx),last(kdx)));
		end
		
              %   quiver(cs(idx).xlim(1)-x0,cs(idx).ylim(1)-y0, ...       
              %                  diff(cs(idx).xlim),diff(cs(idx).ylim),0,'k');
         end 
end

