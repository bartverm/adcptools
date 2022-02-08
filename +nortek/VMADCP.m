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
            nfiles = numel(cur_files);
            
            ad2cp_files = cell(nfiles,1);
            gps_files = cell(nfiles,1);
            for cf = 1:numel(cur_files)
                try
                    dst=unzip(fullfile(cur_files(cf).folder, cur_files(cf).name), unzip_dest);
                catch err
                    error([cur_files(cf).name,' is not a valid SigVM file.'])
                end
                % store gps and ad2cp file names and read them.
                fgps = find(contains(dst, '.anpp','IgnoreCase',true), 1);
                gps_files(cf) = dst(fgps);
                if isempty(fgps)
                    error(['No gps data in file:', cur_files(cf).name])
                end
                fad2cp = find(contains(dst, '.ad2cp','IgnoreCase',true), 1);
                if isempty(fad2cp)
                    error(['No adcp data in file:', cur_files(cf).name])
                end
                ad2cp_files(cf) = dst(fad2cp);
            end
            obj.read_ad2cp(ad2cp_files);
            obj.read_anpp(gps_files);

            % cleanup temporary files
            delete(fullfile(unzip_dest, '*'))
            rmdir(unzip_dest)
        end
        function read_anpp(obj, files)
            obj.gpsraw = obj.read_file(files);
            % anpp processing here
        end
        function val = get_btvel(obj)
        end
        function val = get_bt_vertical_range(obj)	% defined in VMADCP
        end
    end
    properties(Access = protected)
        gpsraw (1, :) uint8
    end
end