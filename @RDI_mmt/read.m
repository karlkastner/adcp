% Fri 14 Sep 18:27:28 CEST 2018
% RDI_mmt :: read
function [first,last, s] = read(filename)
	xml = parseXML(filename);
	s   = xml2struct(xml);
	transect = s.Project.Site_Discharge.Transect;
	%f_C   = fieldnames(s.Project.Site_Discharge.Discharge_Summary.None);
	%k = 0;
	%tag = 'Index_';
	% ensemble indices
	first = [];
	last  = [];
	for idx=1:length(transect)
		first(idx,1) = str2num(transect(idx).Attributes.First); 
		last(idx,1) = str2num(transect(idx).Attributes.Last); 
	%length(f_C)
	%	f = f_C{idx};
	%	if (   length(f)>length(tag) ...
	%	    && strcmp(f(1:length(tag)),tag))
			% ids are zero based, so add 1
			%id = str2num(f(length(tag)+1:end))+1;
			% TODO better parse from b.Project.Site_Discharge.Transect.Attributes.Name
			%eid(id,1) = str2num(s.Project.Site_Discharge.Discharge_Summary.None.(f).StartEnsemble.text.Data);
			%eid(id,2) = str2num(s.Project.Site_Discharge.Discharge_Summary.None.(f).EndEnsemble.text.Data);
	%	end
	%WinRiver.Project.Site_Discharge.Discharge_Summary.None.Index_0.EndEnsemble.text.Data
	%WinRiver.Project.Site_Discharge.Discharge_Summary.BottomTrack.Index_0.StartEnsemble.#text
	%WinRiver.Project.Site_Discharge.Discharge_Summary.GGA.Index_0.StartEnsemble.#text
	%WinRiver.Project.Site_Discharge.Discharge_Summary.VTG.Index_0.StartEnsemble.#text
	%WinRiver.Project.Site_Discharge.Discharge_Summary.GGA2.Index_0.StartEnsemble.#text
	%WinRiver.Project.Site_Discharge.Discharge_Summary.VTG2.Index_0.StartEnsemble.#text
	end
end

