classdef Field < handle & matlab.mixin.Heterogeneous
% Defines a NMEA field in a NMEA message
%
%   Field properties:
%   name - string with field name
%   format - string array with the format of the fields in the NMEA message
%   post_process - post processing function for field data
%
%   Field properties (read only)
%   n_formats - number of format fields to be read
%
%   Field methods (static):
%   default_postproc - default post processing
%   latlong_postproc - latitude, longitude post processing
%   utc_postproc - UTC time post processing
%   mode_postproc - GPS mode post processing
%
%   see also: nmea, nmea.Message
%
%   Copyright 2020, Bart Vermeulen

%     This file is part of the NMEA toolbox.
% 
%     NMEA toolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     NMEA toolbox is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with NMEA toolbox.  If not, see <https://www.gnu.org/licenses/>.

    properties(Constant)
        common_patterns = struct(...
        "float",'(?:\-?\d*\.?\d*)?',... optional minus, one or more digits,0 or 1 dot,0 or more digits 
        "time",'(?:\d{6}\.?\d*)?',... six digits, 0 or 1 dot, 0 or more digits
        "nhex",@(n) ['[a-fA-F0-9]{',num2str(n),'}'],... %any of a to f, A to F, 0 to 9, n times
        "nchr",@(n) ['(?:[a-zA-Z]{',num2str(n),'})?'],... %any of a to z, A to Z, n times
        "ndgt",@(n) ['(?:\-?\d{0,',num2str(n),'})?'],... %optional minus and any of 0-9, n times
        "status",'[AV]?',... % A or V
        "lat",'(?:\d{4}\.?\d*)?,[NS]?',... %four digits, 0 or 1 dot, 0 or mor digits, N or S
        "long",'(?:\d{5}\.?\d*)?,[EW]?') %four digits, 0 or 1 dot, 0 or mor digits, E or W
    end

    properties
        % nmea.Field/name property
        %
        %   scalar string holding name of field in output structure
        %   produced by nmea.Message.parse. Default: "field"
        %
        % see also: nmea, nmea.Field, nmea.Message
        name (1,1) string = "field"
        
        % nmea.Field/format property
        %
        %   scalar or row vector of strings holding fields to be read from
        %   NMEA message. See textscan for how the format of these strings.
        %   Default: "%s"
        %
        %   see also: nmea, nmea.Field, nmea.Message, textscan
        format (1,:) string = "%s"

        pattern (1,:) string = ""
        
        % nmea.Field/post_process property
        %
        %   scalar function handle of function to be called to post-process
        %   data after reading. This function should accept one input which
        %   holds a cell array with one cell element for each format string
        %   given. The function returns the processed data. Default:
        %   @nmea.Field.default_postproc
        %
        %   see also: nmea, nmea.Field, nmea.Field.default_postproc
        post_process (1,1) function_handle = @nmea.Field.default_postproc
    end
    properties (Dependent, SetAccess = private)
        % nmea.Field/n_formats dependent, read only property
        %
        %   n=obj.n_formats()
        %   returns the total number of format fields defined in the Field
        %   objects obj from which the function is called.
        %
        %   see also: nmea, format, nmea.Field
        n_formats

        named_pattern
    end
    methods
        function obj = Field(name, format, pattern, post_process)
            if nargin > 0
                obj.name = name;
            end
            if nargin > 1
                obj.format = format;
            end
            if nargin > 2
                obj.pattern = pattern;
            end
            if nargin > 3
                obj.post_process = post_process;
            end
        end
        function val=get.n_formats(obj)
            val=numel(obj.format);
        end
        function val = get.named_pattern(obj)
            val = obj.get_named_pattern;
        end
    end
    methods(Access = protected)
        function val = get_named_pattern(obj)
            val = "(?<" + obj.name + ">" + obj.pattern + ")";
        end
    end
    methods(Static)
        function out=default_postproc(in)
% returns content of scalar cell, or cell itself when non scalar
%
%   Default post-processing function. If input is a scalar cell, returns
%   the content of the cell. Otherwise returns the cell itself (no change).
%
%   see also: nmea, nmea.Field, latlong_postproc, utc_postproc,
%   mode_postproc
            if isscalar(in)
                out=in{1};
            else
                out=in;
            end
        end

        function out=mode_postproc(in)
% Returns GPSMode array from string input representing the GPS mode
%
%   out = mode_postproc(in) returns an array of type nmea.GPSMode holding
%   the type of GPS fix, given the characters in the input cell in.
%
%   see also: nmea, nmea.Field, latlong_postproc, utc_postproc,
%   mode_postproc, nmea.GPSMode
            out=nmea.GPSMode(ones(size(in{1}),'uint8')*255);
            out(cellfun(@(x) ~isempty(x),in{1}))=nmea.GPSMode(uint8([in{1}{:}]));
        end
    end
end