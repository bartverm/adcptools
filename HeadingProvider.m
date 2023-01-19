classdef HeadingProvider <...
        matlab.mixin.Heterogeneous &...
        matlab.mixin.Copyable & ...
        handle

% Defines interface for classes providing header data
%
%   HeadingProvider methods:
%   has_data - return whether heading data are available
%   heading - return the heading data
%
%   This class allows mixing subclasses in an array. This allows to define
%   a priority of which source of heading is preferred. The first to
%   provide data is used.
%
%   Subclasses must implement the get_has_data and the get_heading methods.
%
%   see also: HeadingProviderInternal, HeadingProviderTFiles
    properties
        heading_misalignment (1,1) double {mustBeFinite, mustBeReal} = 0;

        magnetic_deviation_model (1,1) MagneticDeviationModel = ...
            MagneticDeviationConstant(0);

        % magentic variation
        magnetic_variation (1,1) double {mustBeFinite, mustBeReal} = 0;
    end
    properties(Dependent, SetAccess = private, GetAccess = public)
        magnetic_deviation
    end
    methods (Sealed)
        function val=has_data(obj,adcp)
    % return whether heading data are available
    %
    %   val=has_data(obj,adcp) returns a logical array with the same size 
    %       as obj which is true if corresponding object can return heading
    %       data. adcp must be a scalar ADCP object
    %
    %   see also: HeadingProvider
            if isscalar(obj)
                val=obj.get_has_data(adcp);
            else
                val=false(size(obj));
                for cel=1:numel(obj)
                    val(cel)=obj(cel).get_has_data(adcp);
                end
            end
        end
       
        function val=heading(obj,adcp)
    % returns the heading from the adcp object
    %
    %   val=heading(obj,adcp) return the heading from the provider in the
    %       HeadingProvider array obj. The first object that can provide
    %       heading data will be used.
    %
    % see also: HeadingProvider, ADCP
            if isscalar(obj)
                val=obj.get_heading(adcp) + obj.heading_misalignment;
            else
                idx=find(obj.has_data(adcp),1,"first");
                if isempty(idx)
                    error('No heading data are available')
                end
                val=obj(idx).heading(adcp);
            end
        end
    end
    methods (Abstract, Access=protected)
        val=get_heading(obj,adcp)
        val=get_has_data(obj,adcp)
    end
end