classdef VelocityModel < DataModel

    methods

        function M = get_model(obj, d_time, ~, ~, ~, ~)
            M = ones(numel(d_time), 1, obj.ncomponents);
        end
        function vel = get_velocity(obj, pars, cov_pars,...
                n_bvels, d_time, d_s, d_n, d_z, d_sigma)
            vel = get_data(obj, pars, cov_pars,...
                n_bvels, d_time, d_s, d_n, d_z, d_sigma);
            vel = obj.rotate_matrix(vel, true);
        end
    end
    methods(Access=protected)
        function val = get_ncomponents(~)
            val = 3;
        end
        function val = get_component_names(~)
            val = {'u', 'v', 'w'};
        end
        function val=get_npars(~)
            val = [1 1 1];
        end
        function val=get_names(~)
            val = {{'u'}, {'v'}, {'w'}};
        end
        % function Mrot = rotate_matrix(obj, M, invert)
        %     % assuming M is a n_data x n_pars x n_comp matrix
        %     % Can be done using mat_mult, different implementation here.
        %     % The rotation assumes same sizes of Mu, Mv, Mw.
        %     if nargin < 3
        %         invert = false;
        %     else
        %         assert(islogical(invert) && issscalar(invert))
        %     end
        %     R = obj.rotation_matrix;
        %     if invert
        %         R = R';
        %     end
        %     Mrot = zeros(size(M));
        %     for dim = 1:obj.ncomponents
        %         Mrot(:,:,dim) = R(dim,1)*M(:,:,1) + R(dim,2)*M(:,:,2) + R(dim,3)*M(:,:,3); % Vector quantity
        %     end
        % end

    end
end