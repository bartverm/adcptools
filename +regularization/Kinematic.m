classdef Kinematic < regularization.TaylorBased &...
        regularization.Velocity
    methods(Access = protected)
        function assemble_matrix_private(obj)
            obj.mustBeTaylorModel;
            obj.mustMeetOrderCriteria;
            obj.mustBeVelocityModel;

            % Function that assembles cell-based kinematic boundary conditions.

            zbs = obj.zbsn(1,:);
            zbn = obj.zbsn(2,:);
            dom = obj.domains';

            npars_ne = obj.model.npars_not_expanded;
            assert(all(npars_ne == npars_ne(1)),...
                'Number of non-expanded parameters per component should match')
            npars_ne = npars_ne(1);
            ncells = obj.mesh.ncells;
            npars = sum(obj.model.npars);

            col = [...
                find(obj.find_par(0,'u')),...       u0
                find(obj.find_par(1,'u','sig')),... du/dsig
                find(obj.find_par(0,'v')),...       v0
                find(obj.find_par(1,'v','sig')),... dv/dsig
                find(obj.find_par(0,'w')),...       w0
                find(obj.find_par(1,'w','sig'))];  %dw/dsig
            
            col = reshape(col,npars_ne,ncells,6);

            f_bed = ismember(dom, [6 7 8]);
            n_bed = sum(f_bed);
            f_surf = ismember(dom, [2 3 4]);
            n_surf = sum(f_surf);
            dsig_surf = 1 - obj.mesh.sig_center(f_surf);
            dsig_bed = obj.mesh.sig_center(f_bed);
            dsig_surf = reshape(dsig_surf,1,[]);
            dsig_bed = reshape(dsig_bed,1,[]);

            row = repmat((1 : npars_ne * (n_bed+n_surf))', [1, 6]);
            val = zeros(size(col));

            val(:, f_surf, 5:6) = cat(3,...
                ones(1,n_surf),...
                dsig_surf) .* ones(npars_ne,1);

            val(:, f_bed, :) = repmat(cat(3,...
                zbs(f_bed), ...
               -dsig_bed.*zbs(f_bed), ...
                zbn(f_bed),...
               -dsig_bed.*zbn(f_bed),...
               -ones(1,n_bed),...
                dsig_bed), [size(val,1), 1, 1]);
            
            rhs = zeros(npars_ne, ncells,1);
  
            % add inhomogenous part for surface condition in case of a
            % tidal model
            wl = obj.bathy.water_level;
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
                wl_omega = repmat(wl.model.get_omega,[2 1]);
                wl_omega = wl_omega(:);
                rhs(2:npars_ne,f_surf) = repmat(wl_pars.*wl_omega, [1, n_surf]);
            end
            val = reshape(val,[],6);
            col = reshape(col,[],6);
            rhs = reshape(rhs,[],1);
            is_bad = all(val==0,2);
            val(is_bad, :) = [];
            col(is_bad, :) = [];
            rhs(is_bad) = [];

            obj.C = sparse(...
                row(:),...
                col(:),...
                val(:),...
                size(row,1), ...
                npars*ncells);

            obj.rhs = sparse(rhs);
        end
    end
    methods(Access = protected)
        function val = get_min_order(~)
            val = [...
                0 0 0;... % time
                0 0 0;... % s
                0 0 0;... % n
                0 0 0;... % z
                1 1 1];   % sigma: requires du/dsig, dv/dsig, dw/dsig
        end

        function rhsj = water_level2rhs_element(obj)
            const_names = obj.get_const_names();
            wl = obj.bathy.water_level;
            
            rhsj = zeros([numel(const_names),1]);
            rhsj(1,1) = 0; % subtidal
            if numel(const_names) > 1
                omega = wl.model.get_omega;
                for eq = 1:numel(wl.model.constituents)
                    rhsj(2*eq, 1) = omega(eq)*wl.parameters(2*eq+1);
                    rhsj(2*eq + 1, 1) = -omega(eq)*wl.parameters(2*eq);
                end
            end
        end

    end
end