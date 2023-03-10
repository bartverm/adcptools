function [c, position] = textscan_checked(str,varargin)
	if (isempty(str))
		c = [];
		position = [];
	else
		[c, position] = textscan(str,varargin{:});
	end
end

