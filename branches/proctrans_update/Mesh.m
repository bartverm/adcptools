classdef Mesh < handle
    properties (Dependent)
        ncells
    end
    methods
        function val=get.ncells(obj)
            val=obj.get_ncells();
        end
    end
    methods (Abstract)
        index(obj,n, sigma)
        plot(obj,var)
    end
    methods(Access=protected, Abstract)
        get_ncells(obj)
    end
end