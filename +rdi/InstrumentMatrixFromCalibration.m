classdef InstrumentMatrixFromCalibration < rdi.InstrumentMatrixProvider
% Gets calibrated instrument matrix from pd0 files
%
%   see also: InstrumentMatrixProvider
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            validateattributes(adcp,{'ADCP'},{'scalar'});
            tf=isfield(adcp.raw,'transformation_matrix');
        end
        function i2b=get_i2b_matrix(obj,adcp)
            validateattributes(adcp,{'ADCP'},{'scalar'});
            [b2i, ~, idx]=unique(squeeze(reshape(obj.get_b2i_matrix(adcp),1,adcp.nensembles,[])),'rows'); % get unique matrices
            i2b=zeros(size(b2i)); % initialize output
            for cm=1:size(b2i,1) % for every unique matrix
                i2b(cm,:)=reshape(inv(reshape(b2i(cm,:),4,4)),1,[]);
            end
            i2b=permute(reshape(i2b(idx,:),[],4,4),[4,1,2,3]);
        end
        function b2i=get_b2i_matrix(~,adcp)
            validateattributes(adcp,{'ADCP'},{'scalar'});
            b2i=permute(adcp.int16_to_double(adcp.raw.transformation_matrix)/1e4,[4, 3, 2, 1]);
        end
    end
end