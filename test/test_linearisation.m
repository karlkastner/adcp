
	clf;nf=1;
	% *W'
	plot( [ rmse(meanfilt1(vert(1).model.res_flin,nf)')']);
	ylim([-1.1 1.1]); hold on; set(gca,'ColorOrderIndex',1);
	plot( [ rmse(meanfilt1(vert(1).model.res_f,nf)')'],'--'); 

