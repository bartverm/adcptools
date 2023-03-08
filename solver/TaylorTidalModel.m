classdef TaylorTidalModel < DataModel
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

    %  see also: VelocityModel, TaylorVelocityModel, TidalVelocityModel

    properties
        %         % TaylorExpandedVelocity/s-order order of expansion in s-direction
        %         %
        %         %   1x3 row vector holding non-negative integers that specify up
        %         %   to which order the velocity should be expanded in the
        %         %   s-direction for each of the three cartesian velocity components
        %         %   [2 1 0] means expand x-velocity up to second order in
        %         %   s-direction, y-velocity up to first order and z-velocity up to
        %         %   0th order. For the x-velocity this results in three model
        %         %   parameters (mean x-vel, dx-vel/ds, d^2 x_vel/ds^2), for the
        %         %   y-velocity in two parameters (mean y-vel, dy-vel/dx) and for
        %         %   z-velocity in 1 parameter (mean z-vel).
        %         %
        %         %   Default is: [0, 0, 0]
        %         %
        %         %   see also: TaylorExpandedVelocity, n_order, z_order,
        %         %   sigma_order, time_order
        %         s_order (1,3) double {...
        %             mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
        %
        %         % TaylorExpandedVelocity/n-order order of expansion in n-direction
        %         %
        %         %   1x3 row vector holding non-negative integers that specify up
        %         %   to which order the velocity should be expanded in the
        %         %   n-direction for each of the three cartesian velocity
        %         %   components.
        %         %
        %         %   Default is: [0, 0, 0]
        %         %
        %         %   See for an example: s_order
        %         %
        %         %   see also: TaylorExpandedVelocity, s_order, z_order,
        %         %   sigma_order, time_order
        %         n_order (1,3) double {...
        %             mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
        %
        %         % TaylorExpandedVelocity/z-order order of expansion in z-direction
        %         %
        %         %   1x3 row vector holding non-negative integers that specify up
        %         %   to which order the velocity should be expanded in the
        %         %   z-direction for each of the three cartesian velocity
        %         %   components.
        %         %
        %         %   Default is: [0, 0, 0]
        %         %
        %         %   See for an example: s_order
        %         %
        %         %   see also: TaylorExpandedVelocity, s_order, n_order,
        %         %   sigma_order, time_order
        %         z_order (1,3) double {...
        %             mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
        %
        %         % TaylorExpandedVelocity/sigma-order order of expansion in sigma direction
        %         %
        %         %   1x3 row vector holding non-negative integers that specify up
        %         %   to which order the velocity should be expanded in the
        %         %   sigma-direction for each of the three cartesian velocity
        %         %   components.
        %         %
        %         %   Default is: [0, 0, 0]
        %         %
        %         %   See for an example: s_order
        %         %
        %         %   see also: TaylorExpandedVelocity, s_order, z_order,
        %         %   n_order, time_order
        %         sigma_order (1,3) double {...
        %             mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
        %
        %         % TaylorExpandedVelocity/time-order order of expansion in time
        %         %
        %         %   1x3 row vector holding non-negative integers that specify up
        %         %   to which order the velocity should be expanded in
        %         %   time for each of the three cartesian velocity
        %         %   components.
        %         %
        %         %   Default is: [0, 0, 0]
        %         %
        %         %   See for an example: s_order
        %         %
        %         %   see also: TaylorExpandedVelocity, s_order, z_order,
        %         %   n_order, n_order
        %         %   time_order (1,3) double {...
        %         % mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
        %
        %         % TidalVelocityModel/constituentsU constituents for x-component of velocity
        %         %
        %         %   1xN row vector defining tidal constituents to be included in the model
        %         %   for the x-component of the velocity. Every value given indicates the
        %         %   period of the constituent to be included. For every given value, two
        %         %   model parameters are fitted, which are the coefficients for the sin and
        %         %   cos function. From those amplitude and phases are computed. A residual
        %         %   is always included.
        %         %
        %         %   see also: TidalVelocityModel, constituentsV, constituentsW,
        %         %               get_tidal_pars
        %         constituentsU(1,:) cell = {};
        %
        %         % TidalVelocityModel/constituentsV constituents for y-component of velocity
        %         %
        %         %   1xN row vector defining tidal constituents to be included in the model
        %         %   for the y-component of the velocity. Every value given indicates the
        %         %   period of the constituent to be included. For every given value, two
        %         %   model parameters are fitted, which are the coefficients for the sin and
        %         %   cos function. From those amplitude and phases are computed. A residual
        %         %   is always included.
        %         %
        %         %   see also: TidalVelocityModel, constituentsU, constituentsW,
        %         %               get_tidal_pars
        %         constituentsV(1,:) cell = {};
        %
        %         % TidalVelocityModel/constituentsW constituents for z-component of velocity
        %         %
        %         %   1xN row vector defining tidal constituents to be included in the model
        %         %   for the z-component of the velocity. Every value given indicates the
        %         %   period of the constituent to be included. For every given value, two
        %         %   model parameters are fitted, which are the coefficients for the sin and
        %         %   cos function. From those amplitude and phases are computed. A residual
        %         %   is always included.
        %         %
        %         %   see also: TidalVelocityModel, constituentsU, constituentsV,
        %         %               get_tidal_pars
        %         constituentsW(1,:) cell = {};


        %         TaylorTidalVelocityModel/names names of the different parameter
        %         names. The are coded as

        tidal (1,1) = TidalModel;
        taylor (1,1) = TaylorModel;

    end

    properties(Dependent)
        names
    end

    methods
        function M = get_model(obj, d_t, d_s, d_n, d_z, d_sigma)
            %   [Mu, Mv, Mw]=get_model(obj, time, d_s, d_n, d_z, d_sigma)
            %   procuces the model matrices consisting of Kronecker
            %   products of the Taylor and Tidal velocity models,
            %   respectively.
            T = obj.tidal.get_model(obj,d_t);
            S = obj.taylor.get_model(obj, d_t, d_s, d_n, d_z, d_sigma);

            % Place copies of T as new columns of S: modified Kronecker
            % Size: n_data x n_pars x n_dim
            M = zeros(size(T,1), size(T,2)*size(S,2), size(T,3));

            for dim = 1:obj.ncomponents
                M(:,:,dim) = obj.kron_modified(T(:,:,dim), S(:,:,dim));
            end
        end

        function names = get.names(obj)
            tayl_names = obj.taylor.names;
            tid_names = obj.tidal.names;
            names = cell([obj.ncomponents, 1]);
            for dim = 1:obj.ncomponents
                names{dim,1} = obj.kron_modified_cell(tayl_names{dim, 1}, obj.remove_chars(tid_names{dim, 1}, 1));
            end

        end

    end
    methods(Access=protected)
        function val = get_npars(obj)
            val = ones(1,3) +...
                obj.s_order +...
                obj.n_order +...
                obj.z_order +...
                obj.sigma_order;

            val = val.*[2*numel(obj.constituentsU) + 1, ...
                2*numel(obj.constituentsV) + 1, ...
                2*numel(obj.constituentsW) + 1];
        end
        function val = get_ncomponents(obj)
            val = numel(obj.components);
        end



        function C = stack_matrices(mat_cell)
            if ismatrix(mat_cell)
                C = mat_cell;
            elseif iscell(mat_cell)
                C = zeros(size(mat_cell{1}));
                for i = 1:length(mat_cell)
                    C(:,:,i) = mat_cell{i};
                end
            end
        end

        function C = kron_modified_mat(obj, A, B)
            % Places copies of B in th columns of A after multiplying by
            % those columns, i.e.
            % Can be vectorized but is not the main bottleneck.
            assert((size(A,1)==size(B,1) && ismatrix(A) && ismatrix(B)), 'Design matrix sizes should match')

            C = zeros([size(A,1), size(A,2)*size(B,2)]);
            idx = 1;
            for i = 1:size(A,2)
                col = idx:(idx+size(B,2)-1);
                C(:,col) = repmat(A(:,i), [1, size(B,2)]).*B;
                idx = max(col)+1;
            end
        end

        function C = kron_modified_cell(obj, A, B)
            % Does the same as its matrix equivalent, but now by
            % concatenating the elements of cells.
            assert(((size(A,1)==1) && (size(B,1) == 1)), 'Only row cells allowed')

            % Function is not vectorized, can be done but does not
            % influence program performance
            C = cell([size(A,1), size(A,2)*size(B,2)]);
            idx = 1;
            for i = 1:size(A, 2)
                for j = 1:size(B, 2)
                    C{1,idx} = [A{i}, B{j}];
                    idx = idx + 1;
                end
            end
        end

        function C = remove_chars(obj, C, n)
            % Removes first n characters from every string element
            % belonging to the 1 x N cell C
            % Not vectorized
            for idx = 1:size(C,2)
                C{idx} = C{idx}((n+1):end);
            end
        end
    end

end


