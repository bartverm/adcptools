classdef InstrumentMatrixFromConfString < nortek.InstrumentMatrixProvider
    methods (Access = protected)
        function val = get_has_data(obj, adcp, varargin)
            val = ~isempty(obj.get_xfburst_lines(adcp));
            val = val & isequal([1; 2; 3; 4], unique(adcp.physical_beams_used));
        end
        function val = get_i2b_matrix(obj, adcp, varargin)
            val = obj.get_internal_b2i(adcp);
            val = inv(squeeze(val(1,1,:,:) ) );
            val = shiftdim(val, -2);
            val = repmat(val, [1, adcp.nensembles, 1, 1]);
        end
        function val = get_b2i_matrix(obj, adcp, varargin)
            val = obj.get_internal_b2i(adcp);
            val = cat(3, val(:,:,1:2,:),...
                    cat(4,...
                    val(:,:,3,1)/2, ...
                    val(:,:,4,2)/2, ...
                    val(:,:,3,3)/2, ...
                    val(:,:,4,4)/2),...
                    cat(4,...
                    val(:,:,3,1), ...
                    -val(:,:,4,2), ...
                    val(:,:,3,3), ...
                    -val(:,:,4,4)));
        end  
    end
    methods(Access = protected, Static)
        function val = get_xfburst_lines(adcp)
            conf = adcp.configuration_string;
            lines=strsplit(conf,{'\r','\n'});
            f_lines = ~cellfun(@isempty,regexp(lines,'GETXFBURST'));
            val = lines(f_lines);
        end
        function val = get_internal_b2i(adcp)
            lines = nortek.InstrumentMatrixFromConfString.get_xfburst_lines(adcp);
            if isempty(lines)
                error('No beam to instrument matrix information in configuration string')
            end
            lines = lines{1};
            rows = regexp(lines,'ROWS=(\d)', 'tokens');
            rows = str2double(rows{1}{1});
            cols = regexp(lines,'COLS=(\d)', 'tokens');
            cols = str2double(cols{1}{1});
            val = nan(rows, cols);
            patt = 'M00=(-?\d\.\d+)';
            for cr = 1:rows
                for cc = 1:cols
                    patt(2:3) = [num2str(cr), num2str(cc)];
                    val_string = regexp(lines, patt, 'tokens');
                    val(cr,cc) = str2double(val_string{1}{1});
                end
            end
            val = shiftdim(val, -2);
            val = repmat(val, [1, adcp.nensembles, 1, 1]);
        end
    end
end