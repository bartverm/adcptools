classdef VelocitySolver < ADCPDataSolver
    % Abstract base class to solve ADCP repeat transect velocity data
    %
    %   Subclasses should implement the get_solver_input method.
    %
    %   obj=VelocitySolver(...) allows to set properties upon construction.
    %   Depending on the class objects are assigned to a specific property:
    %   VMADCP -> obj.adcp (required)
    %   Mesh -> obj.mesh (required)
    %   Bathymetry -> obj.bathymetry
    %   XSection -> obj.xs
    %   Filter -> obj.ensemble_filter
    %   VelocityModel -> obj.velocity_model 
    %
    %  VelocitySolver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   ensemble_filter - defines repeat transects to be processed
    %   velocity_model - defines the velocity model to fit
    %
    %   VelocitySolver methods:
    %   get_parameters - get velocity model parameters on mesh
    %   get_velocity - get velocity on mesh
    %   rotate_to_xs - rotates velocity to cross-section direction
    %
    %   see also: VMADCP, Mesh, Bathymetry, XSection, Filter,
    %   TimeBasedVelocitySolver, LocationBasedVelocitySolver



    methods

        function [vel, cov_vel, n_bvels] = get_velocity(obj, pars, cov_pars, n_bvels)
            %   Get velocity from model parameters
            %
            %   vel=get_velocity(obj) computes the velocity in the mesh cells by using
            %   the velocity computed by combining the beam velocities measured at the
            %   same time (this is done by the ADCP class). Velocities that where
            %   measured within a mesh cell (this is determined by the Mesh class) are
            %   averaged.
            %
            %   [vel,cov_vel] = get_velocity(obj) also returns the standard
            %   deviation in the velocity
            if nargin == 1
                [pars, cov_pars, n_bvels] = get_parameters(obj);
            end
            if nargin > 1
                if nargin - 1 < nargout
                    error('VelocitySolver:WrongInputNumber',...
                        'Please pass no input other than object, or as many inputs as required outputs')
                end
            end
            [vel, cov_vel] = ...
                obj.data_model.get_data( ...
                pars, ...
                cov_pars);
        end




        function [vel, cov]=rotate_to_xs(obj, orig_vel, orig_cov)
        % Rotates velocity to direction of cross-section
        %
        %   vel=obj.rotate_to_xs(orig_vel) Rotates the velocity orig_vel
        %   to the direction of the cross-section.
        %
        %   [vel, cov]=obj.rotate_to_xs(orig_vel, orig_cov) also rotates
        %   the covariance matrix
        %
        % See also: VelocitySolver, get_velocity
            [vels, veln] = obj.xs.xy2sn_vel( ...
                orig_vel(:, 1), ...
                orig_vel(:, 2));
            vel = [vels, veln, orig_vel(:, 3)];
            if nargin > 2
                Tsn = obj.xs.xy2sn_tens(orig_cov(:, 1:2, 1:2));
                cov = cat(3, [Tsn, orig_cov(:,3,1:2)],...
                    orig_cov(:,:,3));
            end
        end

    end

end