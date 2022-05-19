classdef TidalVelocityModel < VelocityModel
    properties
        constituentsU(1,:) double {mustBeFinite, mustBePositive} = []
        constituentsV(1,:) double {mustBeFinite, mustBePositive} = []
        constituentsW(1,:) double {mustBeFinite, mustBePositive} = []


    end

    methods
        function [Mu, Mv, Mw] = get_model(obj, d_time, ~, ~, ~, ~) %What about spatial variation?
            % HJ 11-2-22
            % This model fits the following parameters to the velocity
            % within each cell:
            % u = u_0 + sum_n (a_n cos(2pi/T_n * t) + b_n sin(2pi/T_n * t))
            % where n loops over all entered constituents (subtidal
            % constituent is always present)
            % Input:
            % d_time IN HOURS
            % Output:
            % Model matrices such that u = Mu*pars (roughly)
            npars = obj.npars;
            d_secs = seconds(d_time);

            Mu = ones(numel(d_time), npars(1));
            for c = 1:numel(obj.constituentsU)
                Mu(:,2*c) = cos(2*pi/obj.constituentsU(c) * hours(d_secs));
                Mu(:,2*c + 1) = sin(2*pi/obj.constituentsU(c) * hours(d_secs));
            end
            Mv = ones(numel(d_time), npars(2));

            for c = 1:numel(obj.constituentsV)
                Mv(:,2*c) = cos(2*pi/obj.constituentsV(c) * hours(d_secs));
                Mv(:,2*c + 1) = sin(2*pi/obj.constituentsV(c) * hours(d_secs));
            end
            Mw = ones(numel(d_time), npars(3));

            for c = 1:numel(obj.constituentsW)
                Mw(:,2*c) = cos(2*pi/obj.constituentsW(c) * hours(d_secs));
                Mw(:,2*c + 1) = sin(2*pi/obj.constituentsW(c) * hours(d_secs));
            end

        end

        function [pars_h, cov_pars_h] = get_tidal_pars(obj, pars, cov_pars)
            % HJ 11-2-22
            % This function was written to convert estimated harmonic
            % component amplitudes to the standard form:
            % u_n = A_n cos(omega_n t - phi_n) = a_n cos(omega_n t) + b_n
            % sin(omega_n t)
            % Yields A_n = sqrt(a_n^2 + b_n^2)
            % and    phi_n = arctan(b_n / a_n)

            % TO DO: cov_pars_h - analytical formula

            % [pars_h, cov_pars_h] = get_tidal_pars(obj, pars, cov_pars)
            npars = obj.npars;
            subtidal_idx = [1, npars(1) + 1, sum(npars(1:2)) + 1];
            pars_h = zeros(size(pars));
            idx = 1;
            while idx <= size(pars_h,2)
                if ismember(idx, subtidal_idx)
                    pars_h(:,idx) = pars(:,idx);
                    idx = idx + 1;
                else
                    pars_h(:,idx) = sqrt(pars(:,idx).^2 + pars(:,idx+1).^2); % Amplitude of constituent
                    pars_h(:,idx+1) = atan2(pars(:,idx+1), pars(:,idx));     % Phase of constituent
                    %pars_h(:,idx+1) = atan(pars(:,idx+1)./pars(:,idx));

                    idx = idx + 2;
                end
            end

            cov_pars_h = cov_pars; %change this.
        end

        function [vel, cov_vel] = cloud2vel(obj, points, opts)
            % HJ 2-3-22
            % Evaluate velocities on point cloud simultaneously. Could
            % maybe also be incorporated in VelocityModel
            % TO BE WRITTEN
            vel = 0;
            cov_vel = 0;
        end


    end

    methods(Access=protected)
        function val = get_npars(obj)
            val = [2*numel(obj.constituentsU) + 1, ...
                2*numel(obj.constituentsV) + 1, ...
                2*numel(obj.constituentsW) + 1];
        end
    end

end
