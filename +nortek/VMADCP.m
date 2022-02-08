classdef VMADCP < nortek.ADCP & VMADCP
    methods
        function obj = VMADCP(varargin)
            search_string_idx = [];
            search_string = '';
            for ca = 1 : nargin
                if ischar(varargin{ca}) || isStringScalar(varargin{ca})
                    search_string_idx = ca;
                    search_string = varargin{search_string_idx};
                end
            end            
            varargin{search_string_idx} = [];
            obj = obj@nortek.ADCP(varargin{:})
            obj = obj@VMADCP(varargin{:})
            if ~isempty(search_string)
                obj.read_sigvm(search_string)
            end
        end
    end

    methods(Access = protected)
        function read_sigvm(obj, search_string)
            cur_files = dir(search_string);
            if isempty(cur_files)
                error(['Could not find ' ,file_name])
            end
            cur_files([cur_files.isdir]) = [];
            
            unzip_dest = tempname;
            if exist(unzip_dest,"dir")
                error('cannot create temporary directory');
            else
                mkdir(unzip_dest)
            end


            for cf = 1:numel(cur_files)
                try
                    dst=unzip(fullfile(cur_files(cf).folder, cur_files(cf).name), unzip_dest);
                catch err
                    error([cur_files(cf).name,' is not a valid SigVM file.'])
                end
                % store gps and ad2cp file names and read them.
            end
        end
        function val = get_btvel(obj)
        end
        function val = get_bt_vertical_range(obj)	% defined in VMADCP
        end
    end
end