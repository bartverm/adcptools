classdef Data < dynamicprops
    properties (GetAccess=public, SetAccess={?PD0.Reader})
        fixed_leader
        variable_leader
    end
    methods
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