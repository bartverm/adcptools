classdef InstrumentMatrixInternal < InstrumentMatrixProvider
    methods(Access=protected)
        function val = get_has_data(~, adcp, varargin)
            P=inputParser;
            P.addParameter('BottomTracking', false,@(x) isscalar(x) && islogical(x));
            P.parse(varargin{:});

            val = isfield(adcp.raw,'Transformation_Matrices') && ...
                isfield(adcp.raw.Transformation_Matrices,'Frequency') &&...
                isfield(adcp.raw.Transformation_Matrices,'Matrix');
            if P.Results.BottomTracking
                val = val && isfield(adcp.raw,'BottomTrack') &&...
                isfield(adcp.raw.BottomTrack,'BT_Frequency');
            else
                val = val && isfield(adcp.raw,'WaterTrack') &&...
                isfield(adcp.raw.WaterTrack,'WT_Frequency');
            end
        end
        function val = get_i2b_matrix(obj, adcp, varargin)
            val = obj.expand_mats(adcp,...
                pageinv(adcp.raw.Transformation_Matrices.Matrix),...
                varargin{:});
            % unrotated ADCP is downward looking with velocity positive
            % away from the beam
        end
        function val = get_b2i_matrix(obj, adcp, varargin)
            val = obj.expand_mats(adcp,...
                adcp.raw.Transformation_Matrices.Matrix,...
                varargin{:});
        end
        function val = get_beam_orientation_matrix(obj, adcp, varargin)
            val = obj.i2b_matrix(adcp, varargin{:});
            val(:,:,:,4) = [];
            % beam direction is along positive velocity direction,
            % and defined for a downward looking ADCP.
        end
        function val = get_vertical_range_to_cell(~,adcp)
            val = -(adcp.distmidfirstcell+...
                reshape(0:max(adcp.ncells)-1,[],1).*adcp.cellsize);
            % unrotated ADCP is downward looking so vertical range to cell is
            % negative and increases downward.
        end
    end
    methods(Access=protected, Static)
        function val = expand_mats(adcp, mat, varargin)
            % expands the separate matrices for each frequency to the
            % correct matrix in each ensemble, considering the frequency
            % used in that particular ensemble.
            P=inputParser;
            P.addParameter('BottomTracking', false,@(x) isscalar(x) && islogical(x));
            P.parse(varargin{:});

            fid = adcp.raw.file_id;
            all_mat = mat(:,:,:,fid); 
            all_freq = adcp.raw.Transformation_Matrices.Frequency(:,fid); 
            all_freq = shiftdim(all_freq,-2);
            if P.Results.BottomTracking
                freq = shiftdim(adcp.raw.BottomTrack.BT_Frequency, -3);
            else
                freq = shiftdim(adcp.raw.WaterTrack.WT_Frequency, -3);    
            end
            mat_select = repmat(freq == all_freq, [4,4]);
            val = all_mat(mat_select);
            val = reshape(val,4,4,[]);
            val = permute(val,[4,3,1,2]);
        end
    end
end