classdef Sv2SSC < handle
    properties
        adcp (:,1) ADCP
        samples (:,1) acoustics.WaterSample
        averaging_period (1,1) duration = duration(0,3,0)
    end
    methods(Abstract)
        mass_concentration(obj)
        
    end
    methods
        function obj=Sv2SSC(varargin)
            for cv=1:numel(varargin)
                if isa(varargin{cv},'ADCP')
                    obj.adcp=varargin{cv};
                elseif isa(varargin{cv},'acoustics.WaterSample')
                    obj.samples=varargin{cv};
                elseif isa(varargin{cv},'duration')
                    obj.averaging_period=varargin{cv};
                end
            end
        end
        function set.averaging_period(obj,val)
            assert(val>0,'averaging period must be positive')
            obj.averaging_period=val;
        end
        function [adcp_idx, cell_idx, ensemble_idx]=match_adcp_samples(obj)
            n_ensembles=[obj.adcp.nensembles];
            adcp_id=zeros(1,sum(n_ensembles));
            adcp_id(1)=1;
            adcp_id(n_ensembles(1:numel(n_ensembles)-1)+1)=1;
            adcp_id=cumsum(adcp_id);
            dum_n_ensembles=[0 n_ensembles(1:end-1)];
            ens_id=(1:sum(n_ensembles))-dum_n_ensembles(adcp_id);
            sample_time=[obj.samples.time];
            adcp_time=[obj.adcp.time];
            [adcp_idx, cell_idx, ensemble_idx]=deal(cell(size(sample_time)));
            for ct=1:numel(sample_time)
                dt=adcp_time-sample_time(ct);
                if min(abs(dt)) > obj.averaging_period
                    warning(['Sample ',num2str(ct),' has no corresponding ADCP data within given averaging period!'])
                    continue
                end
                f_ens=abs(dt)<obj.averaging_period/2;
                adcp_idx{ct}=adcp_id(f_ens);
                ensemble_idx{ct}=ens_id(f_ens);
            end
            % match cells
            dist=cell(size(obj.adcp));
            for cadcp=1:numel(obj.adcp)
                dist{cadcp}=obj.adcp(cadcp).depth_cell_offset;
            end
            dist=helpers.grow_array(dist{:},'SkipDim',2,'CollectOutput',true);
            dist=cat(2,dist{:});
            dist=-dist(:,:,:,3);
            for cs=1:numel(obj.samples)
                dist_diff=abs(obj.samples(cs).distance-dist(:,ensemble_idx{cs},:));
                [~,cell_idx{cs}]=min(dist_diff,[],1);
            end
        end
    end
end