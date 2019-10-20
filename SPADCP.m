% Wed Dec  3 19:07:32 CET 2014
% Karl Kastner, Berlin
%
%% stream pro acoutic current doppler profiler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
classdef SPADCP < VADCP
	properties
	end % properties
	methods
	function obj = SPADCP(dat,varargin)
		% fix stream pro time
		dat.timeV = dat.timeV1C;
		% stream pro depth pre-filtration
		fdx = find(dat.btrange < 30);
		dat.btrange(fdx) = 0;
		dat.HADCPbeamangle = 20; 

		% recursively call parent constructor
		obj = obj@VADCP(dat,varargin{:},'type',ADCP.SPADCP);

		%obj.time = datenum(dat);
		%dat.timeV1C);
		%obj.type = ADCP.SPADCP;
	end % constructor
	end % methods
end % SPADCP

