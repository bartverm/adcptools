classdef InstrumentMatrixProvider < matlab.mixin.Heterogeneous
% Base class to obtain the instrument calibration matrix
%
%   This class can be subclassed to add a method for providing the
%   Beam to Instrument calibration matrix
%
%   LatLonProvider methods:
%   has_data - has the required data to compute the matrix
%   b2i_matrix - returns the beam 2 instrument matrix
%   i2b_matrix - returns the instrument 2 beam matrix
%
%   see also: VMADCP, InstrumentMatrixFromBAngle,
%   InstrumentMatrixFromCalibration

    methods (Sealed)
        function tf = has_data(obj, adcp)
        % returns whether required data is available to compute matrix
        %
        %   tf = has_data(obj, adcp) returns the array tf which has the
        %   same size as obj, with a logical value indicating whether the
        %   data providers are able to get the instrument matrix from the
        %   given adcp data. adcp must be an object of ADCP class
        %
        %   see also: InstrumentMatrixProvider
            validateattributes(adcp,{'ADCP'},{'scalar'})
            tf=false(size(obj));
            for co=1:numel(obj)
                tf(co)= obj(co).get_has_data(adcp);
            end
        end

        function TM = i2b_matrix(obj, adcp)
        % return instrument to beam calibration matrix
        %
        %   TM = i2b_matrix(obj, adcp) returns the instrument to matrix
        %   calibration matrix as a 1 x nens x 4 x 4 matrix
        %
        %   see also: InstrumentMatrixProvider
            validateattributes(adcp,{'ADCP'},{'scalar'})
            fprovider = find(obj.has_data(adcp), 1);
            if isempty(fprovider)
                error('Cannot compute instrument transformation matrix')
            else
                TM=obj(fprovider).get_i2b_matrix(adcp);
            end
        end
        
        function TM = b2i_matrix(obj, adcp)
        % return the beam to instrument calibration matrix
        %
        %   TM = b2i_matrix(obj, adcp) returns the instrument to matrix
        %   calibration matrix as a 1 x nens x 4 x 4 matrix
        %
        %   see also: InstrumentMatrixProvider
            validateattributes(adcp,{'ADCP'},{'scalar'})
            fprovider = find(obj.has_data(adcp), 1);
            if isempty(fprovider)
                error('Cannot compute instrument transformation matrix')
            else
                TM=obj(fprovider).get_b2i_matrix(adcp);
            end
        end

        function TM = beam_orientation_matrix(obj, adcp)
        % return the beam orientation matrix
        %
        %   TM = beam_orienation_matrix(obj, adcp) returns a matrix
        %   holding the orientation of the beams
        %
        %   see also: InstrumentMatrixProvider
            validateattributes(adcp,{'ADCP'},{'scalar'})
            fprovider = find(obj.has_data(adcp), 1);
            if isempty(fprovider)
                error('Cannot compute instrument transformation matrix')
            else
                TM=obj(fprovider).get_beam_orientation_matrix(adcp);
            end
        end
    end
    methods (Access=protected, Abstract) % methods to implement in subclasses
        get_has_data(adcp)
        get_i2b_matrix(adcp)
        get_b2i_matrix(adcp)
        get_beam_orientation_matrix(adcp)
    end
end