classdef PD0 < handle
    properties
        raw_data
    end
    methods
        function obj=PD0(varargin)
            narginchk(0,1);
            if nargin > 0
                obj.raw_data=varargin{1};
            end
        end
        function time=get_time(obj)
            time=double(obj.raw_data.time_y2k)';
            time(:,2)=time(:,1)*100+time(:,2);
            time(:,1)=[];
            time(:,6)=time(:,6)+time(:,7)/100;
            time(:,7)=[];
            time=datenum(time);
            time=reshape(time,1,[]);
        end
        function ensnum=get_ensnum(obj)
            ensnum=double(obj.raw_data.ensnum_lsb)+double(obj.raw_data.ensnum_msb)*double(intmax('uint16'));
        end
    end
end