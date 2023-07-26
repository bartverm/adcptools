classdef ExternalContinuity < regularization.TaylorBased &...
        regularization.Velocity
    properties
        nscale(1,1) double {mustBeFinite, mustBeReal} = 1;
        sscale(1,1) double {mustBeFinite, mustBeReal} = 1;
    end
    methods (Access = protected)
        function assemble_matrix_private(obj)
            obj.mustBeVelocityModel()
            obj.mustBeTaylorModel()
            obj.mustMeetOrderCriteria()
            assert(isa(obj.mesh, 'SigmaZetaMesh'), ...
                "Only SigmaZetaMesh supported")
            
            wl = obj.bathy.water_level;
            D0 = obj.get_subtidal_depth();
            D0s = -obj.zbsn(1,:);
            D0n = -obj.zbsn(2,:);

            sig_center = reshape(obj.mesh.sig_center,1,[]);
            n_center = obj.mesh.n_middle(obj.mesh.col_to_cell);

            nb = obj.neighbors;
            dom = obj.domains;

            is_in = dom == 0;
            n_in = sum(is_in);
            dsig = sig_center(nb(2, is_in))-sig_center(nb(4, is_in));
            dn = n_center(nb(1, is_in))-n_center(nb(3, is_in));
            sig_center = sig_center(is_in);
            D0 = D0(is_in);
            D0s = D0s(is_in);
            D0n = D0n(is_in);
           
            % number of parameters, when not taylor expanded
            npars_ne = obj.model.npars_not_expanded;
            assert(all(npars_ne == npars_ne(1)),...
                'Number of non-expanded parameters per component should match')
            npars_ne = npars_ne(1);
            ncells = obj.mesh.ncells;
            npars = obj.model.npars;
            npars_total = sum(npars) * ncells;

            par_idx = [...min_order
                find(obj.find_par(1,'u','s')),...
                find(obj.find_par(0,'u')),...
                find(obj.find_par(0,'v')),...
                find(obj.find_par(0,'w'))
                ];
            par_idx = reshape(par_idx, npars_ne, ncells, 4);

            col = cat(3,...
                par_idx(:,       is_in , 1),... du/ds at cell
                par_idx(:, nb(2, is_in), 2),... u top
                par_idx(:, nb(4, is_in), 2),... u bottom
                par_idx(:, nb(1, is_in), 3),... v right
                par_idx(:, nb(3, is_in), 3),... v left
                par_idx(:, nb(2, is_in), 3),... v top
                par_idx(:, nb(4, is_in), 3),... v bottom
                par_idx(:, nb(2, is_in), 4),... w top
                par_idx(:, nb(4, is_in), 4) ... w bottom
                );

            row = repmat((1:npars_ne*n_in)',[1, 9]);
            row = reshape(row,npars_ne,n_in,9);
 
            val = cat(3,...
                 obj.nscale * D0,... du/ds
                 obj.nscale * (1 - sig_center) .* D0s ./ dsig,... du/dsig
                -obj.nscale * (1 - sig_center) .* D0s ./ dsig,... du/dsig
                 obj.sscale * D0 ./ dn,... dv/dn
                -obj.sscale * D0 ./ dn,... dv/dn
                 obj.sscale * (1 - sig_center) .* D0n ./ dsig,... dv/dsig
                -obj.sscale * (1 - sig_center) .* D0n ./ dsig,... dv/dsig
                 obj.sscale * obj.nscale ./ dsig,... dw/dsig
                -obj.sscale * obj.nscale ./ dsig... dw/dsig
                );
            val = repmat(val, [npars_ne, 1, 1]);

            % add correlation with water level variation for tidal model
            % these are added to horizontal derivatives to s and n
            % These are determined in coefficients 1 (du/ds) and 4 & 5
            % (dv/dn)
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
                val(2:npars_ne,:,1) = val(2:npars_ne,:,1) +...
                    obj.nscale*wl_pars; %du/ds
                val(2:npars_ne,:,[4 5]) = val(2:npars_ne,:,[4 5]) * ...
                    (1 + obj.nscale*wl_pars); % dv/dn (different than above
                    % since derivative is computed as a difference)
            end
            obj.C = sparse(...
                row(:), ...
                col(:),...
                val(:),...
                n_in*npars_ne,...
                npars_total);
        end
    end
    methods(Access = protected)
        function val = get_min_order(~)
            val = [...
                0 0 0;... % time
                1 0 0;... % s: du/ds required
                0 0 0;... % n
                0 0 0;... % z
                0 0 0];   % sigma
        end
    end
end