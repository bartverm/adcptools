classdef CoordinateSystem < uint8
% This enumeration defines the four coordinate systems used by ADCPs
    enumeration
        % Coordinates point along the acoustic beam, towards the transducer
        %
        %   This is the most raw coordinate system in which an ADCP
        %   measures velocity. The number of components matches with the
        %   number of beams
        Beam(0)
        
        % Cartesian coodinate system defined with respect to the ADCP
        %
        %   This cartesian coordinate system is defined with respect to the
        %   ADCP and has typically three components (forward, left, up). In
        %   redundant systems (i.e. with more acoustic beams than measured
        %   cartesian component) often an error velocity is added, a
        %   measure for the inhomogeneity of the velocities measured.
        Instrument(1)
        
        % Cartesian coordinate system defined with respect to the Ship
        %
        % This coordinate system is similar to the Instrument coordinate
        % system, but now includes corrections for the tilts of the boat
        % (pitch and roll) and any misalignment of the ADCP and the Ship.
        Ship(2)
        
        % Geographical coordinate system
        %
        % This coordinate system is typically defined with east,north and
        % upward coordinate, and includes also a correction for the heading
        % of the ADCP.
        Earth(3)
    end
end