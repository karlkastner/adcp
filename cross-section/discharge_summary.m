% Fri 28 Jul 10:07:57 CEST 2017
%% compute and store discharge summary
function discharge_summary(opt,base,dt)
%	base = 'sanggau';
%	opt  = sanggau_metadata();

	% hourly data
	% dt   = [];
	% TODO no magic numbers
	% Nlim = 200;

	calib = load_vadcp_discharge(opt);

%	id  = 1:length(opt.filename.vadcp);
%		opt_str = CrossSection.optstr(opt);
%	t = [];
%	Q = [];
%	Q_centre = [];
%	for idx=rvec(id)
%		iname = [opt.filename.vadcp{idx}(1:end-4),opt_str,'.mat'];
%		disp(iname);
%		load(iname);
%
%		t0 = cs.t0;
%		if (isempty(dt))
%			ti = head(t0);
%		else
%			ti = (ceil(t0(1)*24)/24:dt:floor(t0(end)*24)/24)';
%		end
%	
%		N    = cs.N;
%		Q    = [Q; cvec(cs.discharge(ti))];
%
%%		q_tn = cs.q_tn(ti);
%%		ndx       = abs(N)<=Nlim;
%%		Qi        = cs.dw*sum(q_tn)';
%%		Qi_centre = cs.dw*sum(q_tn(ndx,:))';
%
%		t = [t; ti];
%%		Q = [Q; Qi];
%%		Q_centre = [Q_centre; Qi_centre];
%	end % for idx

	t = calib.time;
	Q = calib.Q;

	[t sdx] = sort(t);
	Q = Q(sdx);
	Q_centre = Q;
%	Q_centre = Q_centre(sdx);
	opt_str = CrossSection.optstr(opt);
	obase = ['mat/',base,'-vadcp-summary-',datestr(now,'yyyy-mm-dd-HH-MM'),opt_str];
	save([obase,'.mat'],'t','Q','Q_centre');
	fid = fopen([obase,'.csv'],'w');
	for jdx=1:length(t)
		fprintf(1,'%s; %4g; %4g;\n',datestr(t(jdx),'yyyy/mm/dd HH'), round(Q(jdx)), round(Q_centre(jdx)));
		fprintf(fid,'%s; %4g; %4g;\n',datestr(t(jdx),'yyyy/mm/dd HH'), round(Q(jdx)), round(Q_centre(jdx)));
	end % for jdx
	fclose(fid);

%	fprintf(1,'%s; %4g; %4g;\n',datestr(t.','yyyy/mm/dd HH'), round(Q.'), round(Q_centre.'));
	%opt_str = CrossSection.optstr(opt,'torder',opt.torder);
	%campaign; time; q; q mid; qbs
end % discharge_summary

