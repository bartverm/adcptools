classdef ADCPDataSolver < Solver
    properties
        % Solver/adcp
        %
        %   Scalar VMADCP object holding the adcp data to compute the velocity
        %
        %   see also: Solver, VMADCP
        adcp (1,1) VMADCP

        % Solver/ensemble_filter
        %
        %   Ensemble filter defining ensembles to include
        %
        %   see also: Solver, EnsembleFilter
        ensemble_filter (1,1) EnsembleFilter
    end
    methods
        function obj=ADCPDataSolver(varargin)
            obj = obj@Solver(varargin{:})
            has_vmadcp=false;
            has_bathy=false;
            has_xs=false;
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg,'VMADCP')
                    has_vmadcp=true;
                    obj.adcp=cur_arg;
                elseif isa(cur_arg, 'Bathymetry')
                    has_bathy=true;
                elseif isa(cur_arg,'Filter')
                    obj.ensemble_filter=cur_arg;
                elseif isa(cur_arg,'XSection')
                    has_xs=true;
                end
            end
            if ~has_vmadcp
                error('You must provide a VMADCP object upon construction of a Solver object')
            end
            if ~has_bathy
                obj.bathy=BathymetryScatteredPoints(obj.adcp);
            end
            if ~has_xs
                obj.xs=XSection(obj.adcp);
            end
        end

    end
    methods (Access=protected)
        function [vpos, vdat, xform, time, wl] = get_solver_input(obj)
            vpos = obj.adcp.depth_cell_position; % velocity positions
            vdat = [];
            xform = [];
            time = obj.adcp.time;
            wl = obj.adcp.water_level;

            % get selection of data of current repeat transect and
            % replicate singleton dimensions to get uniform dimensions
            % of inputs
            ens_filt = ~obj.ensemble_filter.bad_ensembles;
            vpos = vpos(:, ens_filt, :,:);
            time = time(ens_filt);
            time = repmat(time, size(vpos, 1), 1, size(vpos, 3));
            wl = wl(ens_filt);
            wl = repmat(wl, size(vpos, 1), 1, size(vpos, 3));
            
            % vectorize
            time = reshape(time, [], 1);
            vpos = reshape(vpos, [], 3);
            wl = reshape(wl, [], 1);
        end
        function [vdat,xform] = filter_and_vectorize(obj,vdat, xform)
            ens_filt = ~obj.ensemble_filter.bad_ensembles;

            % ensemble filter
            xform = xform(:, ens_filt, :, :);
            xform = repmat(xform,...
                 [size(vdat, 1), 1, 1, 1]);
            vdat = vdat(:,ens_filt,:);

            % vectorize
            vdat = reshape(vdat, [], 1);        
            xform = reshape(xform, [], size(xform,4));        
        end
    end

end