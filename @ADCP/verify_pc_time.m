% Thu Mar 12 13:38:37 CET 2015
% Karl Kastner, Berlin
%% verify the time stored in the data file
function [delta, t0, tend] = verify_pc_time(nmea)

	for jdx=find(~cellfun(@isempty,nmea))
		if (isfield(nmea{jdx},'ZDA') ...
			& isfield(nmea{jdx},'RDENS'))
		A = [];
		
		utc  = double(nmea{jdx}.ZDA.UTCtime);
		time = datenum(0,0,0,utc(:,1),utc(:,2),utc(:,3));
		A(nmea{jdx}.ZDA.lineid,1) = time;
		time = double(nmea{jdx}.RDENS.pctime)/Constant.SECONDS_PER_DAY;
		A(nmea{jdx}.RDENS.lineid,2) = time;
		
		id = (1:size(A,1));
		for idx=1:size(A,2)
			fdx = 0 == A(:,idx);
			A(fdx,idx) = interp1(id(~fdx),A(~fdx,idx),id(fdx),'linear','extrap');
		end
		
		delta(jdx) = nanmean(A(fdx,1) - A(fdx,2));
		date       = double(nmea{jdx}.ZDA.date);
		t0(jdx)    = datenum(date(1,1),date(1,2),date(1,3),utc(1,1),utc(1,2),utc(1,3));
		tend(jdx) = datenum(date(end,1),date(end,2),date(end,3),utc(end,1),utc(end,2),utc(end,3));
		%subplot(2,2,1)
		%plot(A)
		%B = [A(:,2).^0 A(:,2).^1];
		%c = B \ A(:,1);
		%A(:,3) = B*c;
		%subplot(2,2,2)
		%plot([A(:,1) A(:,3)])
		end % if
	end % for jdx
end % verify_pc_time

