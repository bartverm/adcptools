classdef ADCP < ADCP
    properties
        % Sontek raw data
        raw(1,1) struct = struct
    end
    methods
        function obj = ADCP(varargin)
            obj = obj@ADCP(varargin{:});
            obj.tilts_provider = sontek.TiltsInternal;
            obj.heading_provider = sontek.HeadingInternal;
            obj.instrument_matrix_provider = ...
                sontek.InstrumentMatrixInternal;
            if nargin>0
                searchstr=varargin{1};
                if ischar(searchstr) ||...
                        (isstring(searchstr) && isscalar(searchstr))
                    files=dir(searchstr); % read all filenames
                    assert(~isempty(files),'readDeployment:NoFileFound',...
                        'Could not find any file'); 
                    files([files.isdir])=[];
                    files=cellfun(@fullfile, {files.folder},...
                        {files.name},'UniformOutput', false);
                elseif iscellstr(searchstr) ||...
                        (isstring(searchstr) && ~isscalar(searchstr))
                    files=searchstr;
                else
                    error('Input should be char or cell cell of char')
                end
                count_good_reads = 0;
                n_files = numel(files);
                tmp_struct = cell(n_files, 1);
                for count_file = 1:n_files
                    if exist(files{count_file}, 'file')
                        try
                            tmp=load(files{count_file});
                        catch err
                            warning(['Could not read file ', ...
                                files{count_file}, ': ',err.message])
                            continue
                        end
                    else
                        warning(['Could not find file: ',...
                            files{count_file}])
                        continue
                    end
                    count_good_reads = count_good_reads + 1;
                    tmp_struct{count_good_reads} = tmp;
                    tmp_struct{count_good_reads}.filename = ...
                        files{count_file};                   
                end
                tmp_struct(count_good_reads + 1 : end)=[];
                obj.raw = obj.cat_structs([tmp_struct{:}]);
                if count_good_reads == 0
                    warning('Could not read any files')
                end
            end
        end
    end
    methods (Access = protected)
        function get_transducer(obj)
        end
        function get_backscatter(obj)
        end
        function get_echo(obj)
        end
        function get_pressure(obj)
        end
        function get_salinity(obj)
        end
        function out = get_temperature(obj)
            out = obj.raw.System.Temperature';
        end
        function out = get_distmidfirstcell(obj)
            out = obj.raw.System.Cell_Start';
        end
        function out = get_time(obj)
            out = datetime(datevec(datenum(2000,1,1)+obj.raw.System.Time/24/3600))';
        end
        function out = get_cellsize(obj)
            out = obj.raw.System.Cell_Size';
        end
        function val = get_beam_angle(obj)
            bom = obj.instrument_matrix_provider.beam_orientation_matrix(obj);
            val = acosd(bom(:,:,:,3));
        end
        function val = get_coordinate_system(obj)
            cor=obj.raw.Setup.coordinateSystem;
            val = repmat(CoordinateSystem.Beam,size(cor));
            val(cor==1) = CoordinateSystem.Ship;
            val(cor==2) = CoordinateSystem.Earth;
        end
        function out = get_ncells(obj)
            out = size(obj.raw.WaterTrack.Velocity,1);
        end
        function out = get_nensembles(obj)
            out = size(obj.raw.System.Time,1);
        end
        function out = get_nbeams(obj)
            out = size(obj.raw.WaterTrack.Velocity,2);
        end
    end
    methods
        function vel = velocity(obj,dst,filter)
        end
        function val = xform(obj,dst,src,varargin)
        end
    end
    methods(Static, Access = private)
        function out = cat_structs(in)
            % concatenates field in different sontek structures
            fields = fieldnames(in);
            for cf = fields'
                current_field = cf{1};
                if isstruct(in(1).(current_field))
                    out.(current_field) = ...
                        sontek.ADCP.cat_structs([in.(current_field)]);
                else
                    if isscalar(in(1).(current_field))
                        out.(current_field) = [in.(current_field)];
                    elseif ischar(in(1).(current_field))
                        out.(current_field) = {in.(current_field)};
                    elseif ismatrix(in(1).(current_field))
                        out.(current_field) = vertcat(in.(current_field));
                    else
                        out.(current_field) = cat(3,in.(current_field));
                    end
                end
            end
            % create a file id to map per-file scalar values to ensembles
            if isfield(in,'System') && all(cellfun(@(x) isfield(x,'Time'),{in.System}))
                n_ens = cellfun(@(x) size(x.Time,1), {in.System});
                out.file_id = zeros(size(out.System.Time,1),1);
                out.file_id(1) = 1;
                out.file_id(cumsum(n_ens(1:end-1))+1) = 1;
                out.file_id = cumsum(out.file_id);
            end
            % fix concatenation of matrices
            if isfield(in,'Transformation_Matrices')
                out.Transformation_Matrices.Frequency = ...
                    reshape(out.Transformation_Matrices.Frequency, 3,[]);
                out.Transformation_Matrices.Matrix = ...
                    reshape(out.Transformation_Matrices.Matrix,4,4,3,[]);
            end
        end
    end
end