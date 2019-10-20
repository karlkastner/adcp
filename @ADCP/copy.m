% Mi 10. Sep 15:26:42 CEST 2014
% Karl Kastner, Berlin
%% copy constructor
function new = copy(obj)
	% Instantiate new object of the same class.
	new = eval(class(obj));
        %new = ADCP(obj.dat);
 
        % Copy all non-hidden properties.
        p = properties(obj);
        for i = 1:length(p)
        	new.(p{i}) = obj.(p{i});
        end
end % copy

