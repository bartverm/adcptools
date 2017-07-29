function [c position] = textscan_checked(str,varargin)
	[c position] = textscan(str,varargin{:});
	if (length(str) ~= position)
		fprintf(1,'Warning: String not parsed until end\n'); 
	end
end

