classdef adcp < handle
    %ADCP Abstract base class defining properties of adcp classes
    %   Details
    
    
    properties(SetAccess={?adcp,?reader})
        n_depth_cells; % Number of depth cells
        n_ensembles; % Number of ensembles
        n_beams; % Number of beams
        beam_matrix; % Transformation matrix for the beams
        orientation; % Pitch roll and heading
        velocity; % velocity data
        echo; % echo data
        correlation; % correlation data
        manufacturer; % additional manufacturer data
    end
    
    
    methods
        
        function obj=adcp(n_cells,n_ensembles,n_beams)
            obj.n_depth_cells=n_cells;
            obj.n_ensembles=n_ensembles;
            obj.n_beams=n_beams;
            
            obj.velocity=nan(obj.n_depth_cells,obj.n_ensembles,obj.n_beams);
            obj.echo=nan(obj.n_depth_cells,obj.n_ensembles,obj.n_beams);
            obj.correlation=nan(obj.n_depth_cells,obj.n_ensembles,obj.n_beams);
        end
        
        function getVelocity(obj,cor,varargin)
            
        end
        function getPosition(obj)
        
        end
    end
    
end

