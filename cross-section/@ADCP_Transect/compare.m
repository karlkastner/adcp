% Sun 16 Sep 16:22:36 CEST 2018
%% discharge summary
function [Q_,tab] = summarise(obj,cmpfile)
	if (nargin() < 2)
	for idx=1:length(obj.time);
		fprintf('%s %4d\r\n', ...
				datestr(obj.time(idx),'yyyy/mm/dd HH:MM'), ...
				round(obj.discharge.centre(idx)));
	end % idx
	end
	
	tab = readtable(cmpfile);
	% eid = tab.EnsembleStart
	% These are default data tags by RDI-WinRiver
	% TODO also for english
	field_C =  {'QGemessen', 'mid';
                    'QOben',     'top';
		    'QUnten',    'bottom';
		    'Flie_geschw_', 'velmag';
		    'Flie_richtung','veldir'};
	for idx=1:size(field_C,1)
		Q.(field_C{idx,2}) = cellfun(@(x) str2num(x),tab.(field_C{idx,1})(2:end-3));
	end
	Q.centre = Q.top + Q.mid + Q.bottom;
	%'mat/2014-04-22-bifurcation-ii-spring/2014-04-22-bifurcation-ii-spring-transect-1.txt')
	% correct for jumpt to next day
	t = cellfun(@(x) datenum(x,'HH:MM:SS'),tab.Startzeit(2:end-3));
	d = cumsum(diff([0;t])<0);
	t = t + d;
	% day
	t = t-floor(t(1)) + floor(obj.time(1));

	% match
	Q_ = struct();
%obj.discharge.centre;
if (1)
	for idx=1:obj.n
		[adt,mdx] = min(abs(t-obj.time(idx)));
		fprintf('%s %4d %s %4d\r\n', ...
				datestr(obj.time(idx),'yyyy/mm/dd HH:MM'), ...
				round(obj.discharge.centre(idx)), ...
				datestr(t(mdx),'yyyy/mm/dd HH:MM'), ...
				round(Q.centre(mdx)) );
		for field_ = {'mid','top','bottom','centre','velmag','veldir'}
			field = field_{1};
			Q_.(field)(idx,1) = Q.(field)(mdx);
		end
		%Q_(idx,2) = Q(mdx);
	end
else
	for idx=1:obj.n
		[adt,mdx] = min(abs(t-obj.time(idx)));
		fprintf('%s %4d %s %4d\r\n', ...
				datestr(obj.time(idx),'yyyy/mm/dd HH:MM'), ...
				round(obj.discharge.centre(idx)), ...
				datestr(t(mdx),'yyyy/mm/dd HH:MM'), ...
				round(Q.centre(mdx)) );
		%Q_(idx,2) = Q(mdx);
	end
end
end

