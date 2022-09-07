classdef HeadTiltInternal < nortek.HeadTiltProvider
    methods(Access = protected)
        function val = get_has_data(~, adcp)
            if ~isa(adcp, 'nortek.ADCP')
                error('HeadTiltInternal only works with nortek.ADCP objects')
            end
            val = any(adcp.burst_has_data(nortek.BurstBit.AHRS));
        end
        function val = get_head_tilt_matrix(~, adcp)
            pitch = adcp.pitch_internal;
            roll = adcp.roll_internal;
            heading = adcp.heading_internal-90;
            zr=zeros(size(heading));
            on=ones(size(heading));
            sh = sind(heading);
            ch = cosd(heading);
            sp = sind(pitch);
            cp = cosd(pitch);
            sr = sind(roll);
            cr = cosd(roll);
            th = cat(3,...
                cat(4,  ch,  sh, zr, zr),...
                cat(4, -sh,  ch, zr, zr),...
                cat(4,  zr,  zr, on, zr),...
                cat(4,  zr,  zr, zr, on));
            tp = cat(3,...
                cat(4,  cp,  -sp.*sr, -cr.*sp, zr),...
                cat(4,  zr,       cr,     -sr, zr),...
                cat(4,  sp,   sr.*cp,  cp.*cr, zr),...
                cat(4,  zr,       zr,      zr, on));
            val = helpers.matmult(th,tp);
        end
    end
end