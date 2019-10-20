% Fri 14 Sep 20:04:32 CEST 2018
% This sets only transect numbers
% TODO this will stop working when rdi change indices
% RDI_mmt :: write
function write(ofilename,pd0pathname,fid,first,last)

	a=what('RDI_mmt');
	templatefilename = [a.path,filesep,'template-raw.mmt'];

	% read root document
	root = xmlread(templatefilename);

	WinRiver = root.item(0).getChildNodes;
	if (~strcmp(WinRiver.getNodeName,'WinRiver')) error('here'); end

	Project=WinRiver.item(1).getChildNodes;
	if (~strcmp(Project.getNodeName,'Project')) error('here'); end

	Site_Discharge=Project.item(5).getChildNodes;
	if (~strcmp(Site_Discharge.getNodeName,'Site_Discharge')) error('here'); end
	
	% the first one is just updated, not copied
	Transect = Site_Discharge.item(1);
	
	for idx=1:length(first)
		if (~strcmp(Transect.getNodeName,'Transect')) error('here'); end
		Transect.setAttribute('First',num2str(first(idx)));
		Transect.setAttribute('Last',num2str(last(idx)));

		File = Transect.getChildNodes.item(1);
		if (~strcmp(File.getNodeName,'File')) error('here'); end
		if (~isempty(fid))
			pd0pathname_ = [pd0pathname(1:end-4),'-',num2str(fid(idx),'%02d'),'.PD0'];
		else
			pd0pathname_ = pd0pathname;
		end
		File.setAttribute('PathName',pd0pathname_);
		File.setAttribute('TransectNmb',num2str(idx-1,'%d'));
		pd0base = basename(pd0pathname_);
		File.getChildNodes.item(0).setNodeValue(pd0base);

		if (idx<length(first))
			Transect = Transect.cloneNode(1);
			Site_Discharge.appendChild(Transect);
			Transect = Site_Discharge.item(Site_Discharge.getLength-1);
		end % if
	end % for

	xmlwrite(ofilename,root);
	% hack to quick fix win-rivers dependence on attribute order
	str = fileread(ofilename);
	str = regexprep(str,'(TransectNmb="[0-9]*")\s*(Type="[0-9]*")','$2 $1');
	%str = strrep(str,'TransectNmb="0" Type="6"','Type="6" TransectNmb="0"');
	%str = strrep(str,'TransectNmb="0" Type="6"','Type="6" TransectNmb="0"');
	fid = fopen(ofilename,'w');
	fprintf(fid,'%c',str);
	fclose(fid);
end % write xml

