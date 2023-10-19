classdef TaylorBased < regularization.Regularization
% Regularization requiring a Taylor model
%
%   This is a abstract class to be subclassed by regularizations that
%   impose prerequisites on estimated derivatives, therefore requiring a
%   taylor expansions based data model.
%
%   regularizaton.TaylorBased properties (read only):
%   min_order - minimum expansion orders required
%   min_time_order - minimum expansion orders required in time
%   min_s_order - minimum expansion order required in s direction
%   min_n_order - minimum expansion order required in n direction
%   min_z_order - minimum expansion order required in z direction
%   min_sig_order - minimum expansion order required in sigma direction
%
%   regularizaton.TaylorBased methods (protected):
%   mustBeTaylorModel - throws error if data model is not Taylor based
%   mustMeetOrderCriteria - throws error if min expansion order are not met
%
%   see also: regularization.Regularization,
%     regularization.InternalContinuity, regularization.ExternalContinuity
    properties(SetAccess=protected, Dependent)
        % regularization.TaylorBased/min_order
        %
        %   minimum required expansion orders to perform regularization.
        %   size: 5 x ncomponents. Rows correspond to time, s, n, z and
        %   sigma coordinates repsectively. Columns correspond to
        %   components of the modelled variable. E.g. suppose we have a
        %   velocity model having three components (u,v,w), the outputs
        %   will be 5 x 3. 
        %   
        %   see also: regularization.TaylorBased
        min_order (5,:) double {mustBeFinite, mustBeReal, mustBeInteger}

        % regularization.TaylorBased/min_time_order
        %
        %   minimum required expansio orders with respect to time
        %   size: 1 x ncomponents. Columns correspond to
        %   components of the modelled variable. E.g. suppose we have a
        %   velocity model having three components (u,v,w), the output
        %   will be 1 x 3. 
        %   
        %   see also: regularization.TaylorBased
        min_time_order (1,:) double {mustBeFinite, mustBeReal,...
            mustBeInteger}

        % regularization.TaylorBased/min_s_order
        %
        %   minimum required expansion orders with respect to the
        %   s-coordinate.
        %   size: 1 x ncomponents. Columns correspond to
        %   components of the modelled variable. E.g. suppose we have a
        %   velocity model having three components (u,v,w), the output
        %   will be 1 x 3. 
        %   
        %   see also: regularization.TaylorBased
        min_s_order (1,:) double {mustBeFinite, mustBeReal, mustBeInteger}

        % regularization.TaylorBased/min_n_order
        %
        %   minimum required expansion orders with respect to the
        %   n-coordinate.
        %   size: 1 x ncomponents. Columns correspond to
        %   components of the modelled variable. E.g. suppose we have a
        %   velocity model having three components (u,v,w), the output
        %   will be 1 x 3. 
        %   
        %   see also: regularization.TaylorBased
        min_n_order (1,:) double {mustBeFinite, mustBeReal, mustBeInteger}

        % regularization.TaylorBased/min_z_order
        %
        %   minimum required expansion orders with respect to the
        %   z-coordinate.
        %   size: 1 x ncomponents. Columns correspond to
        %   components of the modelled variable. E.g. suppose we have a
        %   velocity model having three components (u,v,w), the output
        %   will be 1 x 3. 
        %   
        %   see also: regularization.TaylorBased
        min_z_order (1,:) double {mustBeFinite, mustBeReal, mustBeInteger}

        % regularization.TaylorBased/min_sigma_order
        %
        %   minimum required expansion orders with respect to the
        %   sigma-coordinate.
        %   size: 1 x ncomponents. Columns correspond to
        %   components of the modelled variable. E.g. suppose we have a
        %   velocity model having three components (u,v,w), the output
        %   will be 1 x 3. 
        %   
        %   see also: regularization.TaylorBased
        min_sig_order (1,:) double {mustBeFinite, mustBeReal,...
            mustBeInteger}
        % regularization.TaylorBased/min_n_order
        %
        %   minimum required expansion orders with respect to the
        %   n-coordinate.
        %   size: 1 x ncomponents. Columns correspond to
        %   components of the modelled variable. E.g. suppose we have a
        %   velocity model having three components (u,v,w), the output
        %   will be 1 x 3. 
        %   
        %   see also: regularization.TaylorBased

    end
    methods
        function val = get.min_order(obj)
            val = obj.get_min_order();
        end
        function val = get.min_time_order(obj)
            val = obj.min_order(1,:);
        end
        function val = get.min_s_order(obj)
            val = obj.min_order(2,:);
        end
        function val = get.min_n_order(obj)
            val = obj.min_order(3,:);
        end
        function val = get.min_z_order(obj)
            val = obj.min_order(4,:);
        end
        function val = get.min_sig_order(obj)
            val = obj.min_order(5,:);
        end
        function val = find_par(obj, varargin)
            val = obj.model.find_par(varargin{:});
            val = repmat(val, [obj.mesh.ncells,1]);
        end
    end

    methods(Access = protected)
        function mustBeTaylorModel(obj)
            name = class(obj);
            assert(all(cellfun(@(x) isa(x,"TaylorModel"), {obj.model}),...
                "all"), strcat("Assembling the ", name(16:end), " matrix requires a TaylorModel"));
        end
        function mustMeetOrderCriteria(obj)
            assert(all(obj.model.lumped_orders >= obj.min_order, 'all'),...
                ['Taylor model does not meet minimum required ',...
                 'expansion orders. Check the object''s min_order ',...
                 'properties.'])
        end
    end
    methods(Access = protected, Abstract)
        get_min_order(obj)
    end
end