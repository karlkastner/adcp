% 2016-06-02 18:35:57.104846914 +0200
% Karl Kastner, Berlin
%% plot transects
function plot_transects(cs,adcp, invalid)
		nt = cs.transect.n;
		r = round(sqrt(nt*2/3));
		c = ceil(nt/r);

		if (nargin() < 3)
			invalid = [];
		end

		valid = cs.filter_transects(adcp,invalid);
		nv = sum(valid);

		first = cs.transect.first;
		last  = cs.transect.last;

		UV  = adcp.ens.velocity.earth(:,1:2);
		% rotate to cross section
		dir = cs.dir;
		R =  [dir(1),dir(2);
                     -dir(2),dir(1)];
		UV = UV*R';

		UVb = adcp.btvel.earth(:,1:2);
		UVb = UVb*R';
		H   = adcp.H;
		Q   = cs.discharge;
		A   = cs.area;

		time = adcp.time;

		namedfigure(1,'Discharge');
		clf();
		subplot(2,1,1);
		plot(sort(abs(Q.total(valid))),(1:nv)/(nv+1),'.');
		title('Cumulative distribution');
		subplot(2,1,2);
		%plot(cs.transect.time(valid),[Q.total(valid) Q.left(valid), Q.right(valid)],'.');
		plot(cs.transect.time(valid),[Q.total(valid) Q.centre(valid)+Q.left(valid), Q.centre(valid)+Q.right(valid)],'.');
		hold on
		plot(cs.transect.time(~valid),Q.total(~valid),'ko');
		text(cs.transect.time,double(Q.total),num2str((1:length(valid))'));
		datetick();

		namedfigure(2,'Area');
		clf();
		subplot(2,1,1);
		plot(sort(abs(A.total(valid))),(1:nv)/(nv+1),'.');
		title('Cumulative distribution');
		subplot(2,1,2);
		%plot(cs.transect.time(valid),[Q.total(valid) Q.left(valid), Q.right(valid)],'.');
		plot(cs.transect.time(valid),[A.total(valid) A.centre(valid)+A.left(valid), A.centre(valid)+A.right(valid)],'.');
		text(cs.transect.time,double(A.total),num2str((1:length(valid))'));
		hold on
		plot(cs.transect.time(~valid),A.total(~valid),'ko');
		datetick();

		namedfigure(3,'XY-Coordinates');
		clf();
		plot(adcp.X,adcp.Y,'k.');
		hold on;
		plot(adcp.X(adcp.ens.tid(:,cs.cdx)>0),adcp.Y(adcp.ens.tid(:,cs.cdx)>0),'ko');
		plot(cs.xlim,cs.ylim,'color',[0.5 0.5 0.5],'linewidth',2);
		axis equal

		for idx=1:cs.transect.n
		fdx = (first(idx):last(idx))';
		N   = adcp.N(fdx,cs.cdx);
		T   = adcp.T(fdx,cs.cdx);

		dt  = cvec(cdiff(time(fdx)));

		figure(3);
		plot(adcp.X(fdx),adcp.Y(fdx),'.');

		namedfigure(4,'N-T coordinates');
		subplot(r,c,idx);
		cla
		if (valid(idx))
			plot(N,T,'b.-');
		else
			plot(N,T,'r.-');
		end
		hold on
		plot(N(1),T(1),'bo-');
		xlim([-0.5*cs.dwidth,0.5*cs.dwidth]);
		title(num2str(round([Q.left(idx) Q.centre(idx) Q.right(idx)])))
		axis off
		
		%plot(adcp.N(fdx(1:end-1)),sA*[qf],'.')
		namedfigure(5,'Bottom elevation');
		subplot(r,c,idx);
		cla
		plot(N,-H(fdx),'.');
		xlim([-0.5*cs.dwidth,0.5*cs.dwidth]);
		ylim([min(-H),max(-H)])
%		title(Q)

		namedfigure(6,'Water velocity (Earth reference)');
		subplot(r,c,idx);
		cla
		plot(N,UV(fdx,:),'.');
		xlim([-0.5*cs.dwidth,0.5*cs.dwidth]);

		namedfigure(7,'Boat velocity');
		subplot(r,c,idx);
		plot(N,UVb(fdx,:),'.');
		xlim([-0.5*cs.dwidth,0.5*cs.dwidth]);

		namedfigure(8,'Time step');
		subplot(r,c,idx);
		cla();
		plot(dt,'.');

%		figure(6)
%		plot(adcp.N(fdx(1:end-1)),sA*[qf],'.')
%		hold on

		end
end % plot_transects

