% Sat 15 Sep 14:28:11 CEST 2018
%% export RDI mmt
function obj=export_mmt(obj,mmtname,pd0name,vadcp)
	if (nargin() < 4)
		RDI_mmt.write(mmtname,pd0name,[],obj.first,obj.last);
	else
		[first,first_fid] = vadcp.file_ensemble_index(obj.first);
		[last,last_fid]   = vadcp.file_ensemble_index(obj.last);
		% transects per file
		first_  = [];
		last_   = [];
		%fid_    = [];
		fid_     = []; %first_fid(1);
		fn = vadcp.file_index();
		for idx=1:length(first)
			%if (first_fid(idx) == fid_(end))
				first_(end+1) = first(idx);
				fid_(end+1)   = first_fid(idx);
			%else
				%mmtname_ = [mmtname(1:end-4),'-',num2str(fid-1,'%02d'),'.mmt'];
				%pd0name_ = [pd0name(1:end-4),'-',num2str(fid-1,'%02d'),'.PD0'];
				%RDI_mmt.write(mmtname_,pd0name_,first_-1,last_-1);
			%	first_(end+1) = first(idx);
			%	fid_(end+1)   = first_fid(idx);
			%end
			if (last_fid(idx) == fid_(end))
				last_(end+1) = last(idx);
			else
				last_(end+1) = fn(fid_(end),2)-fn(fid_(end),1);
				%mmtname_ = [mmtname(1:end-4),'-',num2str(fid-1,'%02d'),'.mmt'];
				%pd0name_ = [pd0name(1:end-4),'-',num2str(fid-1,'%02d'),'.PD0'];
				%RDI_mmt.write(mmtname_,pd0name_,first_-1,last_-1);
				% TODO this fails if transect is split in more than 2 files
				first_(end+1) = 1;
				last_(end+1)  = last(idx);
				fid_(end+1)   = last_fid(idx);
			end % if
		end % for
		% last file
		%mmtname_ = [mmtname(1:end-4),'-',num2str(fid-1,'%02d'),'.mmt'];
		%pd0name_ = [pd0name(1:end-4),'-',num2str(fid-1,'%02d'),'.PD0'];
		%RDI_mmt.write(mmtname_,pd0name_,first_-1,last_-1);

		RDI_mmt.write(mmtname,pd0name,fid_-1,first_-1,last_-1);
	end % if
end % export_mmt
