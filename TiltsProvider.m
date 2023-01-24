classdef TiltsProvider < matlab.mixin.Heterogeneous & handle
% Defines interface for classes providing header data
%
%   HeadingProvider methods:
%   has_data - return whether heading data are available
%   pitch - return the pitch data
%   roll - return the roll data
%
%   This class allows mixing subclasses in an array. This allows to define
%   a priority of which source of heading is preferred. The first to
%   provide data is used.
%
%   Subclasses must implement the get_has_data and the get_heading methods.
%
%   see also: TiltsProviderInternal

    methods (Sealed)
        function val=has_data(obj,adcp)
    % return whether pitch and roll data are available
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
       
        function val=pitch(obj,adcp)
    % returns the pitch from the adcp object
    %
    %   val=heading(obj,adcp) return the pitch from the provider in the
    %       TiltsProvider array obj. The first object that can provide
    %       pitch data will be used.
    %
    % see also: TiltsProvider, ADCP
            if isscalar(obj)
                val=obj.get_pitch(adcp);
            else
                idx=find(obj.has_data(adcp),1,"first");
                val=obj(idx).get_pitch(adcp);
            end
        end
        function val=roll(obj,adcp)
    % returns the roll from the adcp object
    %
    %   val=roll(obj,adcp) return the heading from the provider in the
    %       TiltsProvider array obj. The first object that can provide
    %       roll data will be used.
    %
    % see also: TiltsProvider, ADCP
            if isscalar(obj)
                val=obj.get_roll(adcp);
            else
                idx=find(obj.has_data(adcp),1,"first");
                val=obj(idx).get_roll(adcp);
            end
        end
    end
    methods (Abstract, Access=protected)
        val=get_pitch(obj,adcp)
        val=get_roll(obj,adcp)
        val=get_has_data(obj,adcp)
    end
end