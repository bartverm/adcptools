classdef Sv2SSC < handle
% Base class for backscatter to susp sediment concentration calibrations
%
%   obj=acoustics.Sv2SSC() Construct default object
%
%   obj=acoustics.Sv2SSC(objects) Pass objects to assign properties of
%   calibration. Objects of type ADCP are assigned to the 'adcp' property,
%   objects of type acoustics.WaterSample are assigned to the property
%   'samples' and duration objects are assigned to the averaging period
%   object.
%
%   Sv2SSC properties:
%   adcp - ADCP objects used for calibration
%   samples - acoustics.WaterSample objects for calibration
%   averaging_period - duration object to indicate averaging period
%
%   Sv2SSC methods:
%   mass_concentration - mass concentration in g/L (abstract method)
%   match_adcp_samples - match samples with corresponding ADCP data
%   sv_for_calibration - averaged backscatter strength for calibration
%   plot - plot the backscatter strength vs mass concentration
%
%   see also: ADCP, acoustics, Sv2SSC_Power, Sv2SSC_Sassi
    properties
        % acoustics.Sv2SSC/adcp property
        %
        % ADCP objects to be used for calibration
        %
        % see also: Sv2SSC, ADCP
        adcp (:,1) ADCP
        
        % acoustics.Sv2SSC/samples property
        %
        % acoustics.WaterSample objects to be used for calibration
        %
        % see also: Sv2SSC, WaterSample, acoustics
        samples (:,1) acoustics.WaterSample
        
        % acoustics.Sv2SSC/averaging_period property
        %
        % duration object indicating period for backscatter strength
        % averaging. Default is 3 minutes. Duration must be positive.
        %
        % see also: Sv2SSC, duration
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
        function sv_cal=sv_for_calibration(obj)
        % Returns averaged backscatter strength for each water sample
        %
        %   sv_cal=sv_for_calibration(obj) return an array with the same
        %   size as the samples with the average backscatter strength value
        %   for each sample in dB. 
        %
        %   see also: Sv2SSC, ADCP, ADCP/backscatter
            sv=obj.adcp.backscatter;
            [~, cell_idx, ensemble_idx]=obj.match_adcp_samples();
            sv_cal=nan(size(obj.samples));
            for cs=1:numel(obj.samples)
                i1=cell_idx{cs};
                i2=repmat(ensemble_idx{cs},1,1,size(i1,3));
                i3=repmat(shiftdim(1:size(i1,3),-1),1,size(i1,2),1);
                sv_cal(cs)=nanmean(sv(sub2ind(size(sv),i1(:),i2(:),i3(:))));
            end
        end
        function [adcp_idx, cell_idx, ensemble_idx]=match_adcp_samples(obj)
        % Finds ensembles and cells corresponding to water samples
        %
        %   [adcp_idx, cell_idx, ensemble_idx]=match_adcp_samples(obj)
        %   returns indices of adcp objects, depth cells and ensembles for
        %   each of the water samples. The variables are cells with the
        %   same size as the samples property. Each variable holds indices to
        %   the profile data in the adcp objects
        %
        %   see also: Sv2SSC, ADCP, WaterSamples
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
            if numel(dist) > 1
                dist=helpers.grow_array(dist{:},'SkipDim',2,'CollectOutput',true);
            end
            dist=cat(2,dist{:});
            dist=-dist(:,:,:,3);
            for cs=1:numel(obj.samples)
                dist_diff=abs(obj.samples(cs).distance-dist(:,ensemble_idx{cs},:));
                [~,cell_idx{cs}]=min(dist_diff,[],1);
            end
        end
        function plot(obj)
            ssc=[obj.samples.mass_concentration];
            sv=obj.sv_for_calibration;
            plot(sv,ssc,'o')
            set(gca,'yscale','log')
            xlabel('Backscatter strength (dB)')
            ylabel('Sediment concentration (g/L)')
        end
    end
end



