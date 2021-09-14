classdef ADCPVerticalPosition < handle
% Abstract class defining the vertical position of the ADCP
%
%   ADCPVerticalPosition methods:
%   get_vertical_position - returns the ADCP vertical position for the given time
%
% see also: ADCP, ADCPFixedVerticalPosition

    methods (Abstract)
        % ADCPPosition/get_vertical_position
        %
        % Method should be implemented by subclasses and return a 1xN
        % row vector given a scalar ADCP object. The vector holds the 
        % vertical position of the ADCP in m.
        %
        % see also: ADCPHorizontalPosition, ADCP, ADCPFixedVerticalPosition
        pos=get_vertical_position(obj, adcp)
    end
end