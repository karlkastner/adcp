% Sat 15 Sep 14:08:46 CEST 2018
%% ensemble index eid_f with respect to file for ensemble eid
function [eid_f,fid] = file_ensemble_index(obj,eid)
	% file index
	fid  = obj.dat.FileNumber;
	fid  = fid(eid);
	% ensemble index within file
	% index of first ensemble is 1
	file_index = obj.file_index();
	eid_f   = eid - reshape(file_index(fid,1),size(eid)) + 1;
end

