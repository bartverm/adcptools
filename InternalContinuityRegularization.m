classdef InternalContinuityRegularization < TaylorBasedRegularization
% Impose continuity within a mesh cell
    properties
        nscale(1,1) double {mustBeFinite, mustBeReal} = 1;
        sscale(1,1) double {mustBeFinite, mustBeReal} = 1;
    end
    methods(Access=protected)
        function assemble_matrix_private(obj)
            assemble_matrix_private@TaylorBasedRegularization(obj);
            assert(obj.model_is_velocity,...
                "Regularization requires a VelocityModel");
            assert(isa(obj.mesh, 'SigmaZetaMesh'), ...
                "Only SigmaZetaMesh supported")

            min_order = ...
                [0 0 0;... % t
                 1 0 1;... % s
                 0 1 1;... % n
                 0 0 0;... % z
                 0 0 1]; % sigma

            assert(all(obj.model.lumped_orders >= min_order,"all"),...
                ['First order expansion needed for:\n',...
                'u to s and sigma,\n',...
                'v to n and sigma,\n',...
                'w to sigma']);

            % Function that assembles cell-based continuity equation
            wl = obj.bathy.water_level;

            D0 = reshape(obj.get_subtidal_depth(),1,[]);
            % const_names = obj.get_const_names(); % Cell array
            D0s = reshape(-obj.zbsn(1,:),1,[]);
            D0n = reshape(-obj.zbsn(2,:),1,[]);
            sig = reshape(obj.mesh.sig_center,1,[]);

            % number of parameters, when not taylor expanded
            npars_ne = obj.model.npars_not_expanded;
            assert(all(npars_ne == npars_ne(1)),...
                'Number of non-expanded parameters per component should match')
            npars_ne = npars_ne(1);
            ncells = obj.mesh.ncells;
            npars_total = sum(obj.model.npars) * ncells;

            % find indices of parameters in regularization matrix
            col =[find(obj.find_par(1, 'u', 's' )),... % du/ds
                    find(obj.find_par(1, 'u', 'sig')),... % du/dsig
                    find(obj.find_par(1, 'v', 'n'  )),... % dv/dn
                    find(obj.find_par(1, 'v', 'sig')),... % dv/dsig
                    find(obj.find_par(1, 'w', 'sig'))]; % dw/dsig
            col = reshape(col,npars_ne,ncells,5);
            row = repmat((1:npars_ne*ncells)',[1, 5]);
            row = reshape(row,npars_ne,ncells,5);

            % calculate coefficients of parameters in matrix
            val = cat(3,...
                    obj.nscale*D0,... %du/ds
                    obj.nscale*(1-sig).*D0s,... %du/dsig
                    obj.sscale*D0,... % dv/dn
                    obj.sscale*(1-sig).*D0n,... % dv/dsig
                    obj.sscale*obj.nscale*ones(size(sig))); %dw/dsig
            val = repmat(val, [npars_ne, 1, 1]);

            % add correlation with water level variation for tidal model
            % these are added to horizontal derivatives to s and n
            % coefficients 1 and 3 (see above du/ds, dv/dn)
            if isa(obj.model,'TidalModel')
                assert(isa(wl,'VaryingWaterLevel') &&...
                    isa(wl.model,'TidalModel'),...
                    ['A Tidal velocity model also requires a ',...
                    'VaryingWaterLevel with an underlying tidal model'])
                assert(isequal(...
                    wl.model.constituents,...
                    obj.model.constituents),...
                    ['Constituents of water level model ',...
                    'and data model should match'])
                wl_pars = reshape(wl.parameters(2:npars_ne),[],1);
                val(2:npars_ne,:,[1 3]) = val(2:npars_ne,:,[1 3]) +...
                    obj.nscale*wl_pars;
            end

            % build sparse matrix
            obj.C = sparse(...
                row(:),...
                col(:), ...
                val(:), ...
                npars_ne*ncells, ...
                npars_total);

        end
        function val = model_is_velocity(obj)
            val = isa(obj.model,'VelocityModel');
        end
    end
end