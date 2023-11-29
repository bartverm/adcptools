classdef InternalContinuity <...
        regularization.TaylorBased &...
        regularization.Velocity
% Impose continuity (mass conservation) within a mesh cell
%
%   This class uses derivatives of velocity components estimated from data
%   within a mesh cell and weakly imposes continuity. For this
%   regularization to work the following expension requirements must be met:
%   - u must be expanded to at least first order with respect to the 
%     s-coordinate and the sigma-coordinate. 
%   - v must be expanded to at least first order with respect to the 
%     n-coordinate and the sigma-coordinate. 
%   - w must be expanded to at least first order with respect to the 
%     sigma-coordinate. 
%
%   regularization.InternalContinuity properties:
%   nscale - scaling of variables in n-direction
%   sscale - scaling of variables in s-direction
%
%   see also: regularization.Velocity, regularization.ExternalContinuity
    properties
        % regularization.InternalContinuity/nscale
        %
        %   scaling of variables in n-direction. Default value is 1.
        %
        %   see also: ExternalContinuity
        nscale(1,1) double {mustBeFinite, mustBeReal} = 1;

        % regularization.InternalContinuity/sscale
        %
        %   scaling of variables in s-direction. Default value is 1.
        %
        %   see also: ExternalContinuity
        sscale(1,1) double {mustBeFinite, mustBeReal} = 1;
    end
    methods(Access=protected)
        function assemble_matrix_private(obj)
            obj.mustBeVelocityModel()
            obj.mustBeTaylorModel()
            obj.mustMeetOrderCriteria()
            assert(isa(obj.mesh, 'SigmaZetaMesh'), ...
                "Only SigmaZetaMesh supported")

            wl = obj.bathy.water_level;
            D0 = reshape(obj.get_subtidal_depth(),1,[]);
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
            col =[find(obj.find_par(order = 1, component = 'u',...
                    variable = 's' )),... % du/ds
                  find(obj.find_par(order = 1, component = 'u',...
                    variable = 'sig')),... % du/dsig
                  find(obj.find_par(order = 1, component = 'v',...
                    variable = 'n'  )),... % dv/dn
                  find(obj.find_par(order = 1, component = 'v',...
                    variable = 'sig')),... % dv/dsig
                  find(obj.find_par(order = 1, component = 'w',...
                    variable = 'sig'))]; % dw/dsig
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
    
            obj.rhs = sparse([],[],[],npars_ne*ncells,1,0);
        end
    end
    methods(Access = protected)
        function val = get_min_order(~)
            val = ...
                [0 0 0;... % t 
                 1 0 1;... % s: du/ds and du/dsig required
                 0 1 1;... % n: dv/dn and dv/dsig required
                 0 0 0;... % z
                 0 0 1]; % sigma: dw/dsig required
        end
    end
end