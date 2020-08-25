classdef Interpolator < handle
    properties(SetObservable)
        known double {mustBeFinite}=double.empty(3,0)
    end
    properties(Dependent, SetAccess=protected, GetAccess=public)
        n_dims (1,1) double {mustBeFinite,mustBePositive, mustBeInteger}
    end
    methods
        function obj=Interpolator(varargin)
            if nargin > 0
                obj.known=varargin(1);
            end
        end
        function val=get.n_dims(obj)
            val=size(obj.known,1)-1;
        end
    end
    methods (Abstract)
        val=interpolate(obj,query_position)
    end
end