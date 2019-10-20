% Sa 7. Nov 11:35:24 CET 2015
% Karl Kastner, Berlin
%
%% ensemble indices of each file
%
function obj = assign_file(obj)
	FileNumber = obj.dat.FileNumber;
	obj.file(1).n = max(FileNumber);
	for idx=1:obj.file.n
		obj.file.id{idx} = find(idx == obj.dat.FileNumber);
	end
end

