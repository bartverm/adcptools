classdef ADCP < ADCP
    properties
        raw
    end
    methods
        function obj = ADCP(varargin)
            obj = obj@ADCP(varargin{:});
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
                for cf=1:numel(files)
                    if exist(files{cf},'file')
                        try
                            tmp=load(files{cf});
                        catch err
                            warning(['Could not read file ', ...
                                files{cf}, ': ',err.message])
                            continue
                        end
                    else
                        warning(['Could not find file: ', files{cf}])
                        continue
                    end
                    obj(end).raw=tmp;
                    obj(end).raw.filename=files{cf};
                    obj(end+1)=sontek.ADCP; %#ok<AGROW>
                end
                obj(end)=[];
                if isempty(obj)
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
        function get_temperature(obj)
        end
        function get_distmidfirstcell(obj)
        end
        function get_time(obj)
        end
        function get_cellsize(obj)
        end
        function get_beam_angle(obj)
        end
        function get_coordinate_system(obj)
        end
        function get_ncells(obj)
        end
        function get_nensembles(obj)
        end
        function get_nbeams(obj)
        end
    end
    methods
        function velocity(obj)
        end
        function xform(obj)
        end
    end
end