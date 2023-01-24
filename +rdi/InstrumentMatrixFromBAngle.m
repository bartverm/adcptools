classdef InstrumentMatrixFromBAngle < rdi.InstrumentMatrixProvider
% Computes uncalibrated instrument matrix from the beam angle
%
% see also: InstrumentMatrixProvider
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf=all(isfinite(adcp.beam_angle));
        end
        function i2b=get_i2b_matrix(~,adcp)
            bangle=adcp.beam_angle;
            c=adcp.convexity;
            a=sind(bangle);
            b=cosd(bangle);
            d=sqrt(2)*a/2;
            zr=zeros(size(a));
            i2b=cat(3,...
                cat(4,  c.*a,    zr, b,  d),...
                cat(4, -c.*a,    zr, b,  d),...
                cat(4,    zr, -c.*a, b, -d),...
                cat(4,    zr,  c.*a, b, -d));

        end
        function b2i=get_b2i_matrix(~,adcp)
            bangle=adcp.beam_angle;
            c=adcp.convexity;
            a=1./(2*sind(bangle));
            b=1./(4*cosd(bangle));
            d=a./sqrt(2);
            zr=zeros(size(a));
            b2i=cat(3,...
                cat(4, c.*a, -c.*a,    zr,   zr),...
                cat(4,   zr,    zr, -c.*a, c.*a),...
                cat(4,    b,     b,     b,    b),...
                cat(4,    d,     d,    -d,   -d));
        end
    end
end