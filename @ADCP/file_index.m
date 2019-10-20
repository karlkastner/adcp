% Sa 7. Nov 11:35:24 CET 2015
% Karl Kastner, Berlin
%% first and last ensemble index of of a file
function [id, obj] = file_index(obj)
	FileNumber = obj.dat.FileNumber;
	n = max(FileNumber);
	%obj.file(1).n = max(FileNumber);
	id = zeros(n,2);
	for idx=1:n
		id_ = find(idx == obj.dat.FileNumber);
		if (~isempty(id_))
			id(idx,:) = [id_(1),id_(end)];
		end
	end
end

