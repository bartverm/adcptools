classdef InstrumentMatrixFromPS3 < rdi.InstrumentMatrixProvider
% Gets calibrated instrument matrix from ps3 command output
%
%   properties:
%   ps3_matrix - beam to instrument matrix from ps3 command
%
%   see also: InstrumentMatrixProvider
    properties
% InstrumentMatrixFromPS3/ps3_matrix
%
%   Beam to instrument matrix given by the PS3 command
%
%   see also: InstrumentMatrixFromPS3;
        ps3_matrix(4,4) double = eye(4);
        
    end
    methods(Access=protected)
        function tf=get_has_data(obj,~)
            tf=~isequal(obj.ps3_matrix,eye(4));
        end
        function i2b=get_i2b_matrix(obj,adcp)
            validateattributes(adcp,{'ADCP'},{'scalar'});
            i2b=repmat(shiftdim(inv(obj.ps3_matrix),-2),1,adcp.nensembles,1,1);
        end
        function b2i=get_b2i_matrix(obj,adcp)
            validateattributes(adcp,{'ADCP'},{'scalar'});
            b2i=repmat(shiftdim(obj.ps3_matrix,-2),1,adcp.nensembles,1,1);
        end
    end
end