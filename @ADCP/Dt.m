% Thu Jul 17 21:26:14 WIB 2014
% Karl Kastner, Berlin
%
%% projected distance from transducer to cell centres
%% if the instrument is not tilted, this is the vertical distance (depth)
%% between the transducer and cell centres
%% does not account for transducer depth
function Dt = Dt(obj)
	distmidbin1 = obj.distmidbin1;
	nbins       = obj.nbins;
	binsize     = obj.binsize;
	nbins = max(nbins);
	for idx=1:length(binsize)
		Dt(:,idx) = distmidbin1(idx) + (0:nbins-1)'*binsize(idx);
	end
end

