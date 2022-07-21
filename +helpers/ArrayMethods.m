classdef ArrayMethods < handle
% provides helper function to enable array support in methods
    methods (Access = protected)
        function varargout = array_support(obj, func, varargin)
                siz_out = size(obj);
                argout = cell(numel(obj),nargout);
                if ~isempty(varargin)
                    siz_in = cellfun(@size, varargin, ...
                        "UniformOutput",false);
                    assert(isequal(size(obj), siz_in{:}),...
                        'helper.ArrayMethods:NonMatchingInputSize',...
                        'Size of inputs must match size of object array')
                    varargin = cellfun(@(x) reshape(x,[],1), varargin, ...
                        'UniformOutput',false);
                    varargin = [varargin{:}];
                else 
                    varargin = cell.empty(numel(obj),0);
                end
                parfor co = 1:numel(obj)
                    [argout{co,:}] = feval(func, obj(co), varargin{co,:});
                end

                argout = num2cell(argout,1);
                varargout = cellfun(@(x) reshape(x,siz_out), argout, ...
                    'UniformOutput',false);
        end
        function varargout = plot_array(obj,func,varargin)
            siz_out = size(obj);
            hold_stat=get(gca,'nextplot');
            hold on
            argout = cell(numel(obj),nargout);
            for co=1:numel(obj)
                [argout{co,:}]=feval(func,obj(co),varargin{:});
            end
            argout = cellfun(@(x) num2cell(x,1), argout, ...
                'UniformOutput',false);
            varargout = reshape(argout,siz_out);
            set(gca,'NextPlot',hold_stat)
         end
    end
end