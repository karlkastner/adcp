% 2018-05-28 12:01:28.062116619 +0200
%% plot the transect as a line in cartesian coordinates
function plot(obj,xy0,s,varargin)
	% TODO impl. translate and scale
	if (nargin()<2 || isempty(xy0))
		xy0 = [0,0];
	end
	if (nargin()<3 || isempty(s))
		s = 1;
	end
	plot(s*(obj.xlim-xy0(1)),s*(obj.ylim-xy0(2)),'r.-','linewidth',2,varargin{:});
end

