classdef InstrumentMatrixFromBAngle < InstrumentMatrixProvider
% Computes uncalibrated instrument matrix from the beam angle
%
% see also: InstrumentMatrixProvider
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf=all(isfinite(adcp.beam_angle));
        end
        function i2b=get_i2b_matrix(~,adcp)
            bangle=adcp.beam_angle;
            a=sind(bangle);
            b=cosd(bangle);
            zr=zeros(size(a));
            i2b=cat(3,...
                cat(4,  a,    zr,  b,  zr),...
                cat(4, zr,    -a, zr,  b),...
                cat(4, -a,    zr,  b,  zr),...
                cat(4, zr,     a, zr, -b));

        end
        function b2i=get_b2i_matrix(~,adcp)
            bangle=adcp.beam_angle;
            a=1./(2*sind(bangle));
            b=1./(2*cosd(bangle));
            zr=zeros(size(a));
            b2i=cat(3,...
                cat(4,    a,    zr,    -a,   zr),...
                cat(4,   zr,    -a,    zr,    a),...
                cat(4,    b,    zr,     b,    zr),...
                cat(4,   zr,     b,    zr,    b));
        end
        function val = get_beam_orientation_matrix(adcp)

        end
    end
end