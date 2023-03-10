classdef HeadingProvider <...
        matlab.mixin.Heterogeneous &... % Mix subclasses in array
        matlab.mixin.Copyable & ... % Enable copying of object
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
        % offset between heading and actual ADCP heading
        %
        %   heading_misalignment is a scalar double property which
        %   represents the misalignment of the heading with the actual ADCP
        %   heading. This is mainly usefull for external heading devices.
        %
        % see also: HeadingProvider
        heading_misalignment (1,1) double {mustBeFinite, mustBeReal} = 0;

        % model for the local variation in magnetic field 
        % 
        % magnetic_deviation model should be a scalar
        % MagneticDeviationModel object. The default value is a
        % MagneticDeviationConstant model set to a deviation of 0 degrees.
        %
        % see also: MagneticDevaitionModel, MagneticDeviationConstant,
        %   MagneticDeviationTwoCycle
        magnetic_deviation_model (1,1) MagneticDeviationModel = ...
            MagneticDeviationConstant(0);

        % magentic_variation is the offset between the geographic North and
        % the magnetic North
        magnetic_variation (1,1) double {mustBeFinite, mustBeReal} = 0;
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
                val=obj.get_heading(adcp) + obj.heading_corrections(adcp);
            else
                idx=find(obj.has_data(adcp),1,"first");
                if isempty(idx)
                    error('No heading data are available')
                end
                val=obj(idx).heading(adcp);
            end
        end
        function val = magnetic_deviation(obj,adcp)
        % compute magnetic deviation
        %
        %   val = obj.magnetic_deviation(adcp) returns the magnetic
        %   deviation for each ensemble in the adcp object. The deviation
        %   is computed with the model set in the magnetic_deviation_model
        %   property.
        %
        %   see also: HeadingProvider
            val = obj.magnetic_deviation_model.magnetic_deviation(adcp);
        end
        function val = heading_corrections(obj,adcp)
        % compute total heading correction
        %
        %   val = obj.heading_corrections(adcp) returns the total
        %   correction to the heading, i.e. the sum of heading
        %   misalingment, magnetic variation and magnetic deviation for
        %   every ensemble in the adcp object
        %
        %   see also: HeadingProvider
            val =...
                obj.heading_misalignment +...
                obj.magnetic_variation + ...
                obj.magnetic_deviation(adcp);
        end
    end
    methods (Abstract, Access=protected)
        val=get_heading(obj,adcp)
        val=get_has_data(obj,adcp)
    end
end