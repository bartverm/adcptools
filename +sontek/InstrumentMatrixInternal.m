classdef InstrumentMatrixInternal < InstrumentMatrixProvider
    methods(Access=protected)
        function val = get_has_data(~, adcp)
            val = isfield(adcp.raw,'Transformation_Matrices') && ...
                isfield(adcp.raw.Transformation_Matrices,'Frequency') &&...
                isfield(adcp.raw.Transformation_Matrices,'Matrix') &&...
                isfield(adcp.raw,'WaterTrack') &&...
                isfield(adcp.raw.WaterTrack,'WT_Frequency');
        end
        function val = get_i2b_matrix(obj, adcp)
            val = obj.expand_mats(adcp,...
                pageinv(adcp.raw.Transformation_Matrices.Matrix));
        end
        function val = get_b2i_matrix(obj, adcp)
            val = obj.expand_mats(adcp,...
                adcp.raw.Transformation_Matrices.Matrix);
        end
        function val = get_beam_orientation_matrix(obj, adcp)
            val = obj.i2b_matrix(adcp);
            val(:,:,:,4) = [];
            val = -val;
        end
        function val = get_vertical_range_to_cell(~,adcp)
            val = -adcp.distmidfirstcell-...
                reshape(0:max(adcp.ncells)-1,[],1).*adcp.cellsize;
        end
    end
    methods(Access=protected, Static)
        function val = expand_mats(adcp, mat)
            fid = adcp.raw.file_id;
            all_mat = mat(:,:,:,fid); 
            all_freq = adcp.raw.Transformation_Matrices.Frequency(:,fid); 
            all_freq = shiftdim(all_freq,-2);
            wt_freq = shiftdim(adcp.raw.WaterTrack.WT_Frequency, -3);
            mat_select = repmat(wt_freq == all_freq, [4,4]);
            val = all_mat(mat_select);
            val = reshape(val,4,4,[]);
            val = permute(val,[4,3,1,2]);
        end
    end
end