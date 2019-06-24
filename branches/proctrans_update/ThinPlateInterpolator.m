classdef ThinPlateInterpolator < util.Interpolator
   properties
       
   end
   properties(Access=protected)
       tp
   end
   methods
       function obj=ThinPlateInterpolator(varargin)
           [status, errmsg]=license('checkout','curve_fitting_toolbox');
           if ~status
               error(['ThinPlateInterpolator requires the Curve Fitting Toolbox.\n License error:',errmsg]);
           end
           obj=obj@util.Interpolator(varargin{:});
           addlistener(obj,'position','PostSet',@obj.reset_tp)
       end
       function val=interpolate(obj,query_position)
           validateattributes(query_position,{'numeric'},{'2d','real'},'interpolate','query_position',2)
           if isempty(obj.tp)
               disp('constructing spline...')
               obj.tp=tpaps(obj.position',obj.value');
           end
           val=fnval(obj.tp,query_position')';
       end
   end
   methods(Access=protected)
       function reset_tp(obj,varargin)
           obj.tp=[];
       end
   end
   
end