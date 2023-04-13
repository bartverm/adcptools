classdef VelocityModel < DataModel
    properties(Constant)
        components = {'u', 'v', 'w'};
    end
    properties
        % VelocityModel/rotation
        %
        %   Horizontal rotation applied to velocity. Value must be a
        %   scalar, finite and real double.
        %
        % see also: VelocityModel
        rotation(1,1) double {mustBeFinite, mustBeReal} = 0;
    end
    properties(Dependent, SetAccess=private)
        % ncomponents x ncomponents rotation matrix
        rotation_matrix (:,:) double {mustBeFiniate, mustBeReal}
    end
    methods
        function obj = VelocityModel(varargin)
            obj.parse_class_params_inputs(varargin{:})
        end

        function rotation_matrix = get.rotation_matrix(obj)
            rotation_matrix = [cos(obj.rotation), -sin(obj.rotation), 0;...
                sin(obj.rotation), cos(obj.rotation), 0;...
                0, 0, 1];
        end
        function M = get_model(obj, d_time, ~, ~, ~, ~)
            M = ones(numel(d_time), 1, obj.ncomponents);
            M = obj.rotate_matrix(M);
        end
    end
    methods(Access=protected)
        function val=get_ncomponents(~)
            val = 3;
        end
        function val=get_npars(~)
            val = [1 1 1];
        end
        function val=get_names(~)
            val = {'u', 'v', 'w'};
        end
        function Mrot = rotate_matrix(obj, M)
            % assuming M is a n_data x n_pars x n_comp matrix
            % Can be done using mat_mult, different implementation here.
            % The rotation assumes same sizes of Mu, Mv, Mw.

            R = obj.rotation_matrix;
            Mrot = zeros(size(M));
            if numel(R) == 1
                Mrot = M; % Scalar quantity is not rotated
            else
                for dim = 1:obj.ncomponents
                    Mrot(:,:,dim) = R(dim,1)*M(:,:,1) + R(dim,2)*M(:,:,2) + R(dim,3)*M(:,:,3); % Vector quantity
                end
            end
        end

    end
end