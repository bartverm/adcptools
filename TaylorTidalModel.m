classdef TaylorTidalModel < TaylorModel & TidalModel
    %   Velocity model based on Taylor expansions of tidal amplitudes.
    %
    %   @TaylorTidalModel
    %   npars - number of parameters to be estimated
    %
    %   TaylorTidalModel properties:
    %   @TaylorVelocityModel
    %   s_order - expand to given orders in s coordinate
    %   n_order - expand to given orders in n coordinate
    %   z_order - expand to given orders in z coordinate (not recommended)
    %   sigma_order - expand to given orders in sigma coordinate
    %
    %   @TidalModel
    %   constituents
    %
    %   TaylorTidalModel methods:
    %
    %   M = get_model(obj, time, d_s, d_n, d_z, d_sigma)
    %   procuces the model matrices consisting of Kronecker
    %   products of the Taylor and Tidal velocity models,
    %   respectively.
    %
    %   names = get_parameter_names(obj)
    %   produces a cell of strings of length npars containing the names of
    %   the parameters that are estimated. p = [pu pv pw].

    %  see also: VelocityModel, TaylorModel, TidalModel


    methods
        function val = find_par(obj, ord, comp, vars, const)
            all_comps = obj.get_component_names;
            all_vars = obj.var_names;
            all_const = [{'M0'}, obj.constituents];
            if nargin < 2 || isempty(ord)
                ord = [0 1 2];
            end
            if nargin < 3 || isempty(comp)
                comp = all_comps;
            end
            if nargin < 4 || isempty(vars)
                vars = all_vars;
            end
            if nargin < 5 || isempty(const)
                const = all_const;
            end
            f_tayl = find_par@TaylorModel(obj, ord, comp, vars);
            np_tayl = obj.get_npars_tay();
            f_tid = find_par@TidalModel(obj,comp, const);
            val = false(sum(obj.npars),1);
            pos_tayl = 1;
            pos_tayltid = 1;
            pos_tid = 1;
            for cc = 1:obj.ncomponents
                for par = 1:np_tayl(cc)
                    if f_tayl(pos_tayl)
                        val(pos_tayltid+(0:2)) = f_tid(pos_tid+(0:2));
                    end
                    pos_tayltid = pos_tayltid + 3;
                    pos_tayl = pos_tayl + 1; 
                end
                pos_tid = pos_tid + 3;
            end
        end
        function M = get_model(obj, d_t, d_s, d_n, d_z, d_sigma)
            %   [Mu, Mv, Mw]=get_model(obj, time, d_s, d_n, d_z, d_sigma)
            %   procuces the model matrices consisting of Kronecker
            %   products of the Taylor and Tidal velocity models,
            %   respectively.
            T = get_model@TidalModel(obj, d_t);

            % For a TaylorTidal model, dt is irrelevant. Replace by nan

            S = get_model@TaylorModel(obj, d_t, d_s, d_n, d_z, d_sigma);

            % Place copies of T as new columns of S: modified Kronecker
            % Size: n_data x n_pars x n_dim
            M = zeros(size(T,1), size(T,2)*size(S,2), size(T,3));

            for dim = 1:obj.ncomponents
                M(:,:,dim) = helpers.kron_modified_mat(S(:,:,dim), T(:,:,dim));
            end
        end
    end


    methods(Access=protected)
        function names = get_names(obj)
            tayl_names = get_names@TaylorModel(obj);
            tid_names = get_names@TidalModel(obj);
            names = cell([obj.ncomponents, 1]);
            for dim = 1:obj.ncomponents
                names{dim,1} = helpers.kron_modified_cell(tayl_names{dim, 1}, helpers.remove_chars(tid_names{dim, 1}, 1));
            end

        end
        function val = get_npars(obj)
            ntay = get_npars@TaylorModel(obj);
            ntid = get_npars@TidalModel(obj);
            val = ntay.*ntid;
            %             val = ones(1,3) +...
            %                 obj.s_order +...
            %                 obj.n_order +...
            %                 obj.z_order +...
            %                 obj.sigma_order;
            %
            %             val = val.*[2*numel(obj.constituents) + 1, ...
            %                 2*numel(obj.constituents) + 1, ...
            %                 2*numel(obj.constituents) + 1];
        end
        function val = get_ncomponents(obj)
            val = numel(obj.components);
        end
    end
end


