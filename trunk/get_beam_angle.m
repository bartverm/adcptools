function ang=get_beam_angle(adcp)
% retrieve beam angle from adcp structure in degrees

ang=double(adcp.HADCPbeamangle(1));
if ang~=0
    return;
end
switch adcp.sysconf(1,9:10)
    case '00'
        ang=15;
    case '10'
        ang=20;
    case '11'
        ang=30;
    case '01'
        ang=nan;
    otherwise
        error('unhandled beam angle retrieval')
end
end