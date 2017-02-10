classdef ArrayBlock < handle
    properties
        traits=PD0.ArrayBlock_Traits.empty(0,0);
    end
    properties(GetAccess=protected,SetAccess={?PD0.Reader})
        data
    end
    methods
        function disp(obj)
            fprintf('ArrayBlock with %dx%dx%d %s data\n\n',size(obj.data,1),size(obj.data,2),size(obj.data,3),obj.traits.fields.type)
        end
        function B = subsref(obj,S)
            switch S(1).type
                case '.'
                    switch S(1).subs
                        case 'traits'
                            % Reference to A.x and A.y call built-in subsref
                            B = obj.(S.subs);
                        otherwise
                            % Enable dot notation for all properties and methods
                            error('Unknown or private property');
                    end
                case '()'
                    B=obj.data(S.subs{:});                    
            end
        end
        function obj=ArrayBlock(val)
            if nargin > 0
                obj.traits=val;
            end
        end
        function set.traits(obj,val)
            validateattributes(val,{'PD0.ArrayBlock_Traits'},{'scalar'},'set.traits','traits',2)
            obj.traits=val;
        end
    end
end