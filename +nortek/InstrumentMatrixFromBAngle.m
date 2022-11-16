classdef InstrumentMatrixFromBAngle < nortek.InstrumentMatrixProvider
% Computes uncalibrated instrument matrix from the beam angle
%
% see also: InstrumentMatrixProvider
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf=all(isfinite(adcp.beam_angle),'all');
        end
        function i2b=get_i2b_matrix(~,adcp)
            bangle=adcp.beam_angle;
            a=sind(bangle);
            b=cosd(bangle);
            zr=zeros(1,size(a,2));
            i2b=cat(3,...
                cat(4,  a(:,:,1),    zr,  b(:,:,1),  zr),...
                cat(4, zr,    -a(:,:,2), zr,  b(:,:,2)),...
                cat(4, -a(:,:,3),    zr,  b(:,:,3),  zr),...
                cat(4, zr,     a(:,:,4), zr, b(:,:,4)));

        end
        function b2i=get_b2i_matrix(~,adcp)
            bangle=adcp.beam_angle;
            a=1./(2*sind(bangle));
            b=1./(2*cosd(bangle));
            zr=zeros(1,size(a,2));
            b2i=cat(3,...
                cat(4,    a(:,:,1),    zr,    -a(:,:,3),   zr),...
                cat(4,   zr,    -a(:,:,2),    zr,    a(:,:,4)),...
                cat(4,    b(:,:,1),    zr,     b(:,:,3),    zr),...
                cat(4,   zr,     b(:,:,2),    zr,    b(:,:,4)));
        end
    end
end