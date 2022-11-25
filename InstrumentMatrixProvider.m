classdef InstrumentMatrixProvider < matlab.mixin.Heterogeneous
% Base class to obtain instrument specific calibration matrices
%
%   This class can be subclassed to provide manufacturer specific
%   information on ADCP geometry and transformation matrices
%
%   InstrumentMatrixProvider methods:
%   has_data - has the required data to compute the matrix
%   b2i_matrix - returns the beam 2 instrument matrix
%   i2b_matrix - returns the instrument 2 beam matrix
%   beam_orientation_matrix - beam directions for unrotated ADCP
%   vertical_range_to_cells - vertical distance to depth cells
%
%   see also: VMADCP, ADCP, rdi, nortek, sontek

    methods (Sealed)
        function tf = has_data(obj, adcp, varargin)
        % returns whether required data is available to compute matrix
        %
        %   tf = has_data(obj, adcp) returns the array tf which has the
        %   same size as obj, with a logical value indicating whether the
        %   data providers are able to get the instrument matrix from the
        %   given adcp data. adcp must be an object of ADCP class
        %
        %   Options can be specified through 'ParameterName',
        %   ParameterValue specification:
        %   'BottomTrack'
        %   true|{false}
        %   Specifies whether the matrix is neede for Bottom Tracking%
        %
        %   see also: InstrumentMatrixProvider
            validateattributes(adcp,{'ADCP'},{'scalar'})
            tf=false(size(obj));
            for co=1:numel(obj)
                tf(co)= obj(co).get_has_data(adcp, varargin{:});
            end
        end

        function TM = i2b_matrix(obj, adcp, varargin)
        % return instrument to beam calibration matrix
        %
        %   TM = i2b_matrix(obj, adcp) returns the instrument to beam
        %   calibration matrix as a 1 x nens x 4 x 4 matrix
        %
        %   Options can be specified through 'ParameterName',
        %   ParameterValue specification:
        %   'BottomTrack'
        %   true|{false}
        %   Specifies whether the matrix is neede for Bottom Tracking
        %
        %   see also: InstrumentMatrixProvider
            validateattributes(adcp,{'ADCP'},{'scalar'})
            fprovider = find(obj.has_data(adcp, varargin{:}), 1);
            if isempty(fprovider)
                error('Cannot compute instrument transformation matrix')
            else
                TM=obj(fprovider).get_i2b_matrix(adcp, varargin{:});
            end
        end
        
        function TM = b2i_matrix(obj, adcp, varargin)
        % return the beam to instrument calibration matrix
        %
        %   TM = b2i_matrix(obj, adcp) returns the instrument to matrix
        %   calibration matrix as a 1 x nens x 4 x 4 matrix
        %
        %   Options can be specified through 'ParameterName',
        %   ParameterValue specification:
        %   'BottomTrack'
        %   true|{false}
        %   Specifies whether the matrix is neede for Bottom Tracking
        %
        %   see also: InstrumentMatrixProvider
            validateattributes(adcp,{'ADCP'},{'scalar'})
            fprovider = find(obj.has_data(adcp, varargin{:}), 1);
            if isempty(fprovider)
                error('Cannot compute instrument transformation matrix')
            else
                TM=obj(fprovider).get_b2i_matrix(adcp, varargin{:});
            end
        end

        function TM = beam_orientation_matrix(obj, adcp, varargin)
        % return the beam orientation matrix
        %
        %   TM = beam_orientation_matrix(obj, adcp) returns a matrix
        %   holding the orientation of the beams. 
        %   This is a 1 x nens x nbeams x 3 components matrix. The unit
        %   vectors defined in the fourth dimension point away from the
        %   instrument. I.e. if the instrument with zero tilts points
        %   downward, these vectors also point downward.
        %
        %   Options can be specified through 'ParameterName',
        %   ParameterValue specification:
        %   'BottomTrack'
        %   true|{false}
        %   Specifies whether the matrix is neede for Bottom Tracking%
        %
        %   see also: InstrumentMatrixProvider, ADCP
            validateattributes(adcp,{'ADCP'},{'scalar'})
            fprovider = find(obj.has_data(adcp, varargin{:}), 1);
            if isempty(fprovider)
                error('Cannot compute instrument transformation matrix')
            else
                TM=obj(fprovider).get_beam_orientation_matrix(adcp,...
                    varargin{:});
            end
        end
        function vr = vertical_range_to_cell(obj,adcp)
        % the vertical offset from the ADCP to the depth cell
        %
        %   ncells x nensembles x nbeams matrix holding the distance along
        %   the vertical to the depth cells for an unrotated ADCP, with
        %   zero tilts. If the manufacturer defines the ADCP with zero
        %   tilts as upward looking, these values will be positive and
        %   increasing away from the ADCP. If the ADCP with zero tilts is
        %   defined downward looking, the values will be negatively and
        %   decreasing away from the ADCP.
        %
        %   see also: InstrumentMatrixProvider, ADCP
            validateattributes(adcp,{'ADCP'},{'scalar'})
            fprovider = find(obj.has_data(adcp), 1);
            if isempty(fprovider)
                error('Cannot compute instrument transformation matrix')
            else
                vr=obj(fprovider).get_vertical_range_to_cell(adcp);
            end
        end
    end
    methods (Access=protected, Abstract) % methods to implement in subclasses
        get_has_data(obj, adcp, varargin)
        get_i2b_matrix(obj, adcp, varargin)
        get_b2i_matrix(obj, adcp, varargin)
        get_beam_orientation_matrix(obj, adcp, varargin)
        get_vertical_range_to_cell(obj, adcp)
    end
end