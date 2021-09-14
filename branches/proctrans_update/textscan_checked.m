function [c, position] = textscan_checked(str,varargin)
	if (isempty(str))
		c = [];
		position = [];
	else
		[c, position] = textscan(str,varargin{:});
		if (length(str) ~= position)
			fprintf(1,'Warning: String not parsed until end\n'); 
		end
	end
end

